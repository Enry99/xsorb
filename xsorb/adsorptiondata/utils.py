'''
Utility functions for dataclasses serialization and deserialization
'''

def dict_without_none(data):
    '''Return a dictionary excluding keys with None values.'''
    return dict(x for x in data if x[1] is not None)


def convert_generators_to_lists(obj):
    if isinstance(obj, dict): #if dict, call this function recursively for each element
        return {k: convert_generators_to_lists(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)): #if list or tuple, call on each element
        return [convert_generators_to_lists(item) for item in obj]
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes)):
        # generator, excluding strings and bytes
        try:
            return list(obj)
        except: #pylint: disable=bare-except
            return obj
    else:
        return obj
