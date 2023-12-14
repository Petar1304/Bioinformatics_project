
U okviru projekta imate 2 fajla:

cnv.txt, kolone od interesa su chromosome, start, end, cn, cn1, cn2 cn1 i cn2 nisu uvek dostupne,
         cn = cn1 + cn2 izdvojiti regione na kojima je cn=2, ali cn1=1 i cn2=1 (ili nisu dostupni)

.vcf Naći sve varijante unutar VCF fajla koje upadaju u regione koji su izdvojeni iz CNV fajla
     Sračunati VAF (variant allele frequency) za svaku varijantu Odstraniti varijante čija VAF
     vrednost prelazi 60 Korišćenjem jednog od 1D klaster algoritama (HDBScan, KDE, …) grupisati
      varijante u klastere

Izlaz analize treba da budu:
Grafik klasterovanih varijanti
Pozicije centara klastera


