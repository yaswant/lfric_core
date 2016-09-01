{#- This is the skeleton of the configuration feigning module used in unit -#}
{#- tests. The Jinja templating library is used to insert the actual code. -#}
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module feign_config_mod

  use constants_mod, only : {{kinds | sort | join( ', ' )}}

  implicit none

  private
  public :: {{ namelists.keys() | sort | decorate( 'feign_', '_config' ) | join( ', &\n' + ' '*12 ) }}

  integer(i_native), parameter :: temporary_unit = 3

contains

{%- for name, description in namelists | dictsort %}
{%-   set parameters   = description.getParameters() %}
{%-   set enumerations = description.getEnumerations() %}
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{%-   set procedureName = 'feign_' + name + '_config' %}
  subroutine {{procedureName}}( {{arguments[name] | join( ', &\n' + ' '*(15 + procedureName|length) )}} )
{{-'\n'}}
{%- set onlies = ['read_' + name + '_namelist'] %}
{%- if enumerations %}
{%-   for enum, keys in enumerations|dictsort -%}
{%-     do onlies.extend( ['key_from_' + enum, enum + '_from_key'] ) %}
{%-   endfor %}
{%- endif %}
{%- set moduleName = name + '_config_mod' %}
    use {{moduleName}}, only : {{onlies | join( ', &\n' + ' '*(17 + moduleName|length) )}}

    implicit none
{{-'\n'}}
{%- for param in arguments[name] %}
{%-   set fortranType = parameters[param] %}
{%-   if fortranType.typex == 'character' %}
    character(*), intent(in) :: {{param}}
{%-   else %}
    {{fortranType.typex}}({{fortranType.kind}}), intent(in) :: {{param}}
{%-   endif %}
{%- endfor %}

    integer(i_native) :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_{{name}}_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&{{name}}")' )
{%- for param in arguments[name] %}
{%-   set fortranType = parameters[param] %}
{%-   if param in enumerations %}
    write( temporary_unit, '("{{param}} = ''", A, "''")' ) key_from_{{param}}( {{param}} )
{%-   else %}
{%-     if fortranType.typex=='logical' %}
{%-       set formatString = '", L' %}
{%-     elif fortranType.typex=='integer' %}
{%-       set formatString = '", I0' %}
{%-     elif fortranType.typex=='real' %}
{%-       set formatString = '", E14.7' %}
{%-     elif fortranType.typex=='character' %}
{%-       set formatString = '\'\'", A, "\'\'"' %}
{%-     endif %}
    write( temporary_unit, '("{{param}} = {{formatString}})' ) {{param}}
{%-   endif %}
{%- endfor %}
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_{{name}}_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'feign_{{name}}_config: Unable to close temporary file'

  end subroutine feign_{{name}}_config
{%- endfor %}

end module feign_config_mod
