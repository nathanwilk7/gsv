def add(a,b):
    return a+b

def sub(a,b):
    return a-b

funcs = []
funcs.append(add)
funcs.append(sub)
for func in funcs:
    print(func(4,3))
