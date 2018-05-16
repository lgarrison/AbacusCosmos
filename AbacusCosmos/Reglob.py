from glob import glob
import os
import os.path
import re

def reglob(pattern, regex, return_capture=False):
    """
    A replacement for glob that allows regex filtering.
    Instances of '{}' in the pattern will first be replaced by '*'
    to get all possible matches, then the matches will be filtered
    by replacing '{}' with `regex`.
    
    Parameters
    ----------
    pattern: str
        The globbing pattern
    regex: str, or tuple of str, or dict of str
        The regex(es)
    return_capture: int, list of int, or bool, optional
        Return the nth capture group matched by the regex.
        The list of matches will instead be a list of captures.
        Default: False.
        
    Returns
    -------
    matches: list
        The sorted list of matches
    """
    pattern = os.path.normpath(pattern)

    try:
        len(return_capture)
        flat_capture = False
    except:
        flat_capture = True
        return_capture = [return_capture]
    return_capture = [int(r) for r in return_capture if r not in (False, None)]
    
    if type(regex) is str:
        regex = (regex,)
    nreg = len(regex)
    
    if type(regex) is tuple:
        glob_pat = pattern.format(*(('*',)*nreg))  # all {} are replaced with *
        re_pat = pattern.format(*regex)
    elif type(regex) is dict:
        glob_pat = pattern.format(**{k:'*' for k in regex})  # all {key:} are replaced with *
        re_pat = pattern.format(**regex)
    res = sorted(f for f in glob(glob_pat) if re.match(re_pat, f))
    if return_capture:
        for i,r in enumerate(res):
            res[i] = [re.match(re_pat, r).group(rc) for rc in return_capture]
            if flat_capture:
                res[i] = res[i][0]
    return res
