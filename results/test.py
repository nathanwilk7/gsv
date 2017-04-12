class Test:
    def __init__(self):
        self.a = 0
        self.b = 0

def modif(test):
    test.a = 1

test = Test()
#print(test._replace(b=15))
modif(test)
print(test.a)
#print(test)
