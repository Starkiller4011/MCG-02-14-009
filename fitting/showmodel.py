#!/usr/bin/env python3
# Helper functions
def showmodel(m):
    '''
    Print current model information
    '''
    print('Model: {}'.format(m.expression))
    print()
    print('{:4} {:4} {:12} {:10} {:7} {:15} {:12}'.format('P#',
                                                          'C#',
                                                          'Component',
                                                          'Parameter',
                                                          'Unit',
                                                          'Value',
                                                          'Errors'))
    print('--'*38)
    pid = 1
    for cid, component in enumerate(m.componentNames):
        for parameter in eval('m.{}.parameterNames'.format(component)):
            u = eval('m.{}.{}.unit'.format(component, parameter))
            val = eval('m.{}.{}.values[0]'.format(component, parameter))
            err = eval('m.{}.{}.error[:2]'.format(component, parameter))
            print('{:<4} {:<4} {:<12} {:<10} {:<7} {:<10.5} \
                   ({:<10.5}, {:<10.5})'.format(pid, cid + 1, component,
                                                parameter, u, val,
                                                err[0], err[1]))
            pid += 1
