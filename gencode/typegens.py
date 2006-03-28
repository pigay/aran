import os.path
import re

_pat="\@%s\@"

def scantypedefs(file):
    g={}
    l={}
    execfile(file,g,l)
    if not l.has_key('types_list'):
        raise "types_list not found in file \"%s\"." % file

    tmp = l['types_list']
    res = []
    for type in tmp:
        d = {}
        for k,v in type.items():
            d[re.compile(_pat % k)] = v
        res.append(d)
    return res

def scansequence(s,repl_dict):
    res=s
    for pat,repl in repl_dict.items():
        res = pat.sub(repl,res)
    return res

def scanfile(file,repl_dict, dir):
    dstfile=scansequence(file,repl_dict)
    sf=open(file,'r')
    df=open(os.path.join (dir, os.path.basename (dstfile)),'w')
    for l in sf.readlines():
        df.write(scansequence(l,repl_dict))

    sf.close()
    df.close()


if __name__ == "__main__":
    import sys

    arglist = sys.argv[1:]

    templates = []
    files = []
    dir = os.getcwd ()
    for a in arglist:
        if a.startswith ("-f="):
            templates.append(a[3:])
        elif a.startswith ("-dir="):
            dir = a[5:]
        else:
            files.append(a)

    types = []

    for t in templates:
        types.extend(scantypedefs(t))

    for t in types:
        for f in files:
            scanfile(f,t, dir)
