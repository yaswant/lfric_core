##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

from transformations import Dynamo0p3ColourTrans, DynamoOMPParallelLoopTrans

def trans(psy):
    ctrans = Dynamo0p3ColourTrans()
    otrans = DynamoOMPParallelLoopTrans()

    # Loop over all of the Invokes in the PSy object
    for invoke in psy.invokes.invoke_list:

        print "Transforming invoke '{0}' ...".format(invoke.name)
        schedule = invoke.schedule

        # Colour loops unless they are on W3
        for loop in schedule.loops():
            if loop.field_space != "w3":
                schedule, _ = ctrans.apply(loop)

        # Add OpenMP to loops unless they are over colours
        for loop in schedule.loops():
            if loop.loop_type != "colours":
                schedule, _ = otrans.apply(loop)

        # take a look at what we've done
        schedule.view()

    return psy
