from array import *
from si import prime_factors
rKlasse = 0
neutral = 0
elemente = 0
eMents = []
eMentsCount = []


"""
Module zur Verwendung von Gruppen, immer mit start_gruppen aufrufen
"""
def start_gruppen():
    global rKlasse, neutral, elemente
    rKlasse = int(input("Restklasse = "))
    neutral = int(input("Neutrales Element = "))
    elemente = rKlasse
    gruppe = int(input("1 für Addidative | 2 für Multiplikative :"))
    lookFor = int(input("Suche nach Ordnung von Element?(y=1|no=0) "))
    if (lookFor == 1):
        finish = int(input("How many are you looking for?: "))
        for l in range(finish):
            pruef(gruppe, lookFor)
    else:
        another = int(input("Check alle Elemente + Ordnung einer Gruppe (1) | Alle Elementordnungen + deren Anzahl (2) "))
        pruef(gruppe, lookFor, another)


def pruef(gruppe, lookFor, another):
    if (another == 2):
        if (gruppe == 1):
            if (lookFor == 1):
                numb = int(input("Search for = "))
                searchADD(numb)
            elif(lookFor == 0):
                addidative()
            else:
                print ("Invalid Search Input")
                start_gruppen()

        elif (gruppe == 2):
            if (lookFor == 1):
                numb = int(input("Search for = "))
                searchMULT(numb)
            elif(lookFor == 0):
                multiplikative()
            else:
                print ("Invalid Search Input")
                start_gruppen()
        else:
            print ("Invalid Group")
            print ("Input Valid One")
            start_gruppen()
    elif (another == 1):
        if(gruppe == 2 and lookFor == 0):
            multi_all()
        elif(gruppe == 2 and lookFor == 1):
            return False
        elif (gruppe == 1 and lookFor == 0):
            add_all()
        elif (gruppe == 1 and lookFor == 1):
            return False
            
    else:
        print("Invalid Checker")
        start_gruppen()


def searchADD(numb):
    tmp = 0
    ergebnis = 0
    for x in range(elemente):
        tmp += numb
        tmp %= rKlasse
        eMents.append(tmp)
        if (tmp == neutral):
            ergebnis = (x+1)
            break
        else:
            continue
    print ("Elemente:")
    for t in range(len(eMents)):
        formated = ("[%d] | %d" % (t+1, eMents[t]))
        print (formated)
    print ("Ordnung: ", ergebnis)
    return 0

def searchMULT(numb):
    tmp = 0
    for x in range(1,elemente):
        tmp = pow(numb,x,rKlasse)
        eMents.append(tmp)
        if(tmp == neutral):
            ergebnis = x
            break
        else:
            continue
    print ("Ordnung: ", ergebnis)
    print ("Elemente: ")
    for t in range(len(eMents)):
        formated = ("[%d] | %d " % (t+1, eMents[t]))
        print (formated)
    return 0



def addidative():

    for i in range(elemente):
        tmp = 0
        for x in range(elemente):
            tmp += i
            tmp %= rKlasse
            if (tmp == neutral):
                if (eMents.count(x) == 0):   
                    eMents.append(x)
                    eMentsCount.append(x)
                else:
                    eMentsCount.append(x)
                    
                break
            else:
                continue
    
    for k in range(len(eMents)):
        formated = ("Ordnung: %d |  Vorkommnis dieser Ordnung: %d" % (eMents[k]+1, eMentsCount.count(eMents[k])))
        print (formated)

def multiplikative():
    for i in range(1,elemente):
        for x in range(1,elemente):
            tmp = pow(i,x, rKlasse)
            if (tmp == neutral):
                if (eMents.count(x) == 0):   
                    eMents.append(x)
                    eMentsCount.append(x)
                else:
                    eMentsCount.append(x)
                    
                break
            else:
                continue
    
    for k in range(len(eMents)):
        formated = ("Ordnung: %d |  Vorkommnis dieser Ordnung: %d" % (eMents[k], eMentsCount.count(eMents[k])))
        print (formated)



def multi_all():
    for i in range(1,elemente):
        list_elem = []
        for x in range(1,elemente):
            tmp = pow(i,x, rKlasse)
            list_elem.append(tmp)
            if (tmp == neutral):
                if (eMents.count(x) == 0):   
                    eMents.append(x)
                    eMentsCount.append(x)
                else:
                    eMentsCount.append(x)
                String = ("Ordnung: %d" % x)
                print(list_elem, String)
                break
            else:
                continue
    
    for k in range(len(eMents)):
        formated = ("Ordnung: %d |  Vorkommnis dieser Ordnung: %d" % (eMents[k], eMentsCount.count(eMents[k])))
        print (formated)

def add_all():
    for i in range(elemente):
        list_elem = []
        tmp = 0
        for x in range(elemente):
            tmp += i
            tmp %= rKlasse
            list_elem.append(tmp)
            if (tmp == neutral):
                if (eMents.count(x) == 0):   
                    eMents.append(x)
                    eMentsCount.append(x)
                else:
                    eMentsCount.append(x)
                String = ("Ordnung: %d" % (x+1))
                print(list_elem, String)
                break
            else:
                continue
    
    for k in range(len(eMents)):
        formated = ("Ordnung: %d |  Vorkommnis dieser Ordnung: %d" % (eMents[k]+1, eMentsCount.count(eMents[k])))
        print (formated)

"""
def ordnung_check(g, n, neutral):
    if (pow(g, n)):
        for i in prime_factors(n, n):
            if(pow(g, (n/i))==neutral):
                return False
    else:
        return False
    return True
"""
