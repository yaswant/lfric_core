{#- This is the skeleton of the namelist loading module.                   -#}
{#- The Jinja templating library is used to insert the actual code.        -#}
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the {{listname}} namelist.
!>
module {{listname}}_config_mod

  use constants_mod, only : {{kindlist | join( ', ' )}}
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use ESMF,          only : ESMF_VM, ESMF_VMBroadcast, ESMF_SUCCESS

  implicit none

  private
  public ::
{%- if enumerations %}
{%-   for enumeration in enumerations.keys() | sort %}
{{- ' ' }}{{-enumeration}}_from_key, key_from_{{enumeration}}{{', &\n           '}}
{%-   endfor %}
{%- endif -%}
{{' '}}read_{{listname}}_namelist, {{listname}}_is_loadable, {{listname}}_is_loaded

{%- if enumerations %}
{{-'\n'}}
{%-   for enumeration, pairs in enumerations | dictsort %}
{%-     for pair in pairs %}
  integer(i_native), public, parameter :: {{listname}}_{{enumeration}}_{{pair.key}} = {{pair.value}}
{%-     endfor %}
{%-   endfor %}
{%- endif %}

{%- if parameters %}
{{-'\n'}}
{%-   for name, ftype in parameters | dictsort %}
  {{ftype.typex}}({{ftype.kind}}), public, protected :: {{name}}
{%-   endfor %}
{%- endif %}

  logical :: namelist_loaded = .false.

{%- if enumerations %}
{{-'\n'}}
{%-   for enumeration, pairs in enumerations | dictsort %}
  character(str_def), parameter :: {{enumeration}}_key({{pairs | length()}}) &
{%         set indent = '          = [character(len=str_def) :: ' %}
{%-     for pair in pairs %}
{%-       if not loop.first %}
{%-         set indent = ' ' * indent | length() -%}
, &{{'\n'}}
{%-        endif %}
{{- indent }}'{{ pair.key }}'
{%-     endfor %}]
{%-   endfor %}
{%- endif %}

contains

{%- for enumeration, pairs in enumerations | dictsort %}

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> \param[in] key Enumeration key.
  !>
  integer(i_native) function {{enumeration}}_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    key_index = 1
    do
      if (trim({{enumeration}}_key(key_index)) == trim(key)) then
        {{enumeration}}_from_key = key_index + {{listname}}_{{enumeration}}_{{pairs[0]['key']}} - 1
        return
      else
        key_index = key_index + 1
        if (key_index > ubound({{enumeration}}_key, 1)) then
          write( log_scratch_space, &
                 '("Key ''", A, "'' not recognised for {{listname}} {{enumeration}}")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function {{enumeration}}_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> \param[in] value Enumeration value.
  !>
  character(str_def) function key_from_{{enumeration}}( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: key_index

    key_index = value - {{listname}}_{{enumeration}}_{{pairs[0]['key']}} + 1
    if (key_index < lbound({{enumeration}}_key, 1) &
        .or. key_index > ubound({{enumeration}}_key, 1)) then
      write( log_scratch_space, &
             '("Value ", I0, " is not in {{listname}} {{enumeration}}")' ) value
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_from_{{enumeration}} = {{enumeration}}_key( key_index )

  end function key_from_{{enumeration}}
    {%- endfor %}

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \param [in] file_unit Unit number of the file to read from.
  !> \param [in] vm ESMF VM object of current run.
  !> \param [in] local_rank Rank of current ESMF process.
  !>
  subroutine read_{{listname}}_namelist( file_unit, vm, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    type(ESMF_VM),     intent(in) :: vm
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, vm, local_rank
    {%- if enumerations -%}
    , {{enumerations.keys() | sort | join( ', ' )}}
    {%- endif -%}
    {{' '}})
  end subroutine read_{{listname}}_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, vm, local_rank
    {%- if enumerations -%}
    , {{ enumerations.keys() | sort | decorate( 'dummy_' ) | join( ', ' ) }}
    {%- endif -%}
    {{' '}})

    use constants_mod, only : RMDI, IMDI{%- if constants %}, {{constants | join( ', ' )}}{%- endif %}

    implicit none

    integer(i_native), intent(in) :: file_unit
    type(ESMF_VM),     intent(in) :: vm
    integer(i_native), intent(in) :: local_rank
    {%- if enumerations %}
    {%-   for enumeration, pairs in enumerations.iteritems() %}
    integer(i_native), intent(out) :: dummy_{{enumeration}}
    {%-   endfor %}
    {{-'\n'}}
    {%-   for enumeration in enumerations.keys() | sort %}
    character(str_def) :: {{enumeration}}
    {%-   endfor %}
    {%- endif %}

    namelist /{{listname}}/ {{ variables.keys() | sort | join( ', &\n' + ' '*(16+listname|length) ) }}

    integer(i_native) :: condition

    {%- for kind in kindlist | sort %}
      {%- if kindcounts[kind] > 0 %}
        {%- if kind[0] == "i" %}
    integer({{kind}}) :: 
        {%- elif kind[0] == "r" %}
    real({{kind}}) :: 
        {%- elif kind[0] == "l" %}
    integer(i_native) :: 
        {%- elif kind[0] == "s" %}
    character({{kind}}) :: 
        {%- endif 
         %} bcast_{{kind}}({{kindcounts[kind]}})
      {%- endif %}
    {%- endfor %}

    {%-for name, ftype in parameters | dictsort %}
    {%-  if ftype.kind[0] == "i" and name not in enumerations %}
    {{name}} = IMDI
    {%-  elif ftype.kind[0] == "r" %}
    {{name}} = RMDI
    {%-  elif ftype.kind[0] == "l" %}
    {{name}} = .FALSE.
    {%-  elif ftype.kind[0] == "s" or name in enumerations %}
    {{name}} = ""
    {%-  endif %}
    {%- endfor %}

    if (local_rank == 0) then

      read( file_unit, nml={{listname}}, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      {{-'\n'}}      
      {%- for name in enumerations.keys() | sort %}
      dummy_{{name}} = {{name}}_from_key( {{name}} )
      {%- endfor %}
      {%- for kind in kindlist | sort %}
      {%-   set count = 1 %}
      {%-   for name, ftype in parameters | dictsort %}
      {%-     if ftype.kind == kind %}
      bcast_{{kind}}({{count}})
      {%-       if name in enumerations 
                 %} = dummy_{{name}}
      {%-       elif name in logicals
                 %} = merge(0, 1, {{name}})
      {%-       else
                 %} = {{name}}
      {%-       endif %}
      {%-     set count = count + 1 %}
      {%-     endif %}
      {%-   endfor %}
      {%- endfor %}

    end if
    {{-'\n'}}
    {%- for kind in kindlist | sort %}

    {%-   if kindcounts[kind] > 0 %}
    call ESMF_VMBroadcast( vm, bcast_{{kind}}, {{kindcounts[kind]}} 
      {%- if kind[0] == "s" %}*{{kind}}{%- endif %}, 0, rc=condition )
    if (condition /= ESMF_SUCCESS) then
      write(log_scratch_space, "(A)") &
          "Failed to broadcast {{kindcounts[kind]}} {{kind}} varaible/s"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    {%-   endif %}
    {%- endfor %}
    
    if (local_rank /= 0) then
      {{-'\n'}}      
      {%- for kind in kindlist | sort %}
      {%-   set count = 1 %}
      {%-   for name, ftype in parameters | dictsort %}
      {%-     if ftype.kind == kind %}
      {%-       if name in enumerations %}
      dummy_{{name}} = 
      {%-       elif name in logicals %}
      {{name}} = (
      {%-       else %}
      {{name}} =
      {%-       endif 
                     %} bcast_{{kind}}({{count}}) {%
                if name in logicals %}== 0 ) {%
                endif %}
      {%-     set count = count + 1 %}
      {%-     endif %}
      {%-   endfor %}
      {%- endfor %}

    end if

    namelist_loaded = .true.
{% for name, code in initialisation.iteritems() %}
    {{name}} = {{code[0]}}
{% endfor %}
  end subroutine read_namelist

  !> Can this namelist be loaded?
  !>
  !> \return True if it is possible to load the namelist.
  !>
  function {{listname}}_is_loadable()

    implicit none

    logical :: {{listname}}_is_loadable

    {{listname}}_is_loadable = .not. namelist_loaded

  end function {{listname}}_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \return True if the namelist has been loaded.
  !>
  function {{listname}}_is_loaded()

    implicit none

    logical :: {{listname}}_is_loaded

    {{listname}}_is_loaded = namelist_loaded

  end function {{listname}}_is_loaded

end module {{listname}}_config_mod
