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
    return_capture: int, optional
        Return the nth capture group matched by the regex.
        The list of matches will instead be a list of captures.
        Default: False.
        
    Returns
    -------
    matches: list
        The sorted list of matches
    """
    # First remove repeated '/', which will not be returned by glob
    # and thus will fail the regex match
    #pattern = re.sub(r'(/)\1+', r'\1', pattern)
    pattern = os.path.normpath(pattern)
    
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
        res = [re.match(re_pat, r).group(return_capture) for r in res]
    return res