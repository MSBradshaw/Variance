
def split_file(name):
    file = open('pileups-by-species/' + name,'r')
    outits1 = open('its1/' + name,'w')
    outits2 = open('its2/' + name,'w')
    out58s = open('5.8s/' + name,'w')
    count = 0
    for line in file:
        if count > 600 and count < 792:
            outits1.write(line)
        elif count >= 792 and count  < 951:
            out58s.write(line)
        else:
            outits2.write(line)
        count = count + 1
    outits1.close()
    outits2.close()
    out58s.close()

split_file('R.hayd-78.pileup')

names = ['R.hayd-35.pileup',  'R.mel-10.pileup',    'R.poly-3O.pileup',  'R.port-8A.pileup',    'R.shus-2.pileup','R.hayd-78.pileup',  'R.mel-68_B.pileup',  'R.poly-3T.pileup',  'R.port-95.pileup',    'R.shus-3.pileup','R.mel-00.pileup',   'R.pari-03.pileup',   'R.poly-5B.pileup',  'R.port-96.pileup',    'R.shus-4.pileup','R.mel-01.pileup',   'R.pari-05.pileup',   'R.poly-8G.pileup',  'R.port-96TA.pileup',  'R.shus-5.pileup','R.mel-02.pileup',   'R.pari-65.pileup',   'R.port-3A.pileup',  'R.port-97.pileup',    'R.shus-6.pileup']
for name in names:
    split_file(name)