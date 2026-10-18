/* The genes the drawings show. Every one is real: a window of AMAP run on a
   genome whose right answer is known (tests/benchmark), with every Bakta call
   in it ("bakta", with the case the run gave it), every placement of a
   reference protein the classification saw ("ref", "piece" marking the other
   piece of a gene a contig end cuts), Bakta's known false positives aligned
   there ("bad"), how the case 1 calls' proteins match the
   reference genome ("refMatch"), what the classification wrote ("corrected")
   and the genes AMAP kept ("amap"). Where check_cds.py changed a gene, the
   window carries the stretch of genome it read ("translate"); a window with a
   "kind" or "exact" carries what that step's check needs. tools/check_site.py
   runs every window through bin/ and fails if anything comes out different.
   Written by tools/figure_data.py.

   Everything after the "=" is JSON. */
window.AMAP_FIGURE_DATA = {
 "windows": {
  "hero": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 1,
   "to": 12400,
   "bakta": [
    {
     "tag": "L4_CDC1551_00001",
     "s": 34,
     "e": 1524,
     "st": "+",
     "name": "dnaA",
     "product": "Chromosomal replication initiator protein DnaA",
     "case": 4
    },
    {
     "tag": "L4_CDC1551_00002",
     "s": 2052,
     "e": 3260,
     "st": "+",
     "name": "dnaN",
     "product": "DNA polymerase III (beta chain) DnaN (DNA nucleotidyltransferase)",
     "case": 2
    },
    {
     "tag": "L4_CDC1551_00003",
     "s": 3280,
     "e": 4437,
     "st": "+",
     "name": "recF",
     "product": "DNA replication and repair protein RecF (single-strand DNA binding protein)",
     "case": 2
    },
    {
     "tag": "L4_CDC1551_00004",
     "s": 4482,
     "e": 4997,
     "st": "+",
     "name": "",
     "product": "Conserved hypothetical protein",
     "case": 4
    },
    {
     "tag": "L4_CDC1551_00005",
     "s": 5123,
     "e": 7267,
     "st": "+",
     "name": "gyrB",
     "product": "DNA gyrase (subunit B) GyrB (DNA topoisomerase (ATP-hydrolysing)) (DNA topoisomerase II) (type II DNA topoisomerase)",
     "case": 4
    },
    {
     "tag": "L4_CDC1551_00006",
     "s": 7302,
     "e": 9818,
     "st": "+",
     "name": "gyrA",
     "product": "DNA gyrase (subunit A) GyrA (DNA topoisomerase (ATP-hydrolysing)) (DNA topoisomerase II) (type II DNA topoisomerase)",
     "case": 2
    },
    {
     "tag": "L4_CDC1551_00007",
     "s": 9963,
     "e": 10160,
     "st": "+",
     "name": "",
     "product": "hypothetical protein",
     "case": 4
    },
    {
     "tag": "L4_CDC1551_00008",
     "s": 10523,
     "e": 10828,
     "st": "+",
     "name": "",
     "product": "hypothetical protein",
     "case": 4
    },
    {
     "tag": "L4_CDC1551_00009",
     "s": 11123,
     "e": 11341,
     "st": "-",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    },
    {
     "tag": "L4_CDC1551_00010",
     "s": 11421,
     "e": 11528,
     "st": "-",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    },
    {
     "tag": "L4_CDC1551_00011",
     "s": 11555,
     "e": 11692,
     "st": "-",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    },
    {
     "tag": "L4_CDC1551_00012",
     "s": 11874,
     "e": 12311,
     "st": "-",
     "name": "",
     "product": "putative membrane protein",
     "case": 2
    }
   ],
   "ref": [
    {
     "prot": "Rv0001",
     "name": "dnaA",
     "s": 1,
     "e": 1524,
     "st": "+",
     "len": 508,
     "match": 1524,
     "num": 1524,
     "complete": true
    },
    {
     "prot": "Rv0002",
     "name": "dnaN",
     "s": 2052,
     "e": 3260,
     "st": "+",
     "len": 403,
     "match": 1209,
     "num": 1209,
     "complete": true
    },
    {
     "prot": "Rv0003",
     "name": "recF",
     "s": 3280,
     "e": 4437,
     "st": "+",
     "len": 386,
     "match": 1155,
     "num": 1158,
     "complete": true
    },
    {
     "prot": "Rv0004",
     "name": "Rv0004",
     "s": 4434,
     "e": 4997,
     "st": "+",
     "len": 188,
     "match": 564,
     "num": 564,
     "complete": true
    },
    {
     "prot": "Rv0005",
     "name": "gyrB",
     "s": 5240,
     "e": 7267,
     "st": "+",
     "len": 676,
     "match": 2028,
     "num": 2028,
     "complete": true
    },
    {
     "prot": "Rv0006",
     "name": "gyrA",
     "s": 7302,
     "e": 9818,
     "st": "+",
     "len": 839,
     "match": 2508,
     "num": 2517,
     "complete": true
    },
    {
     "prot": "Rv0007",
     "name": "Rv0007",
     "s": 9914,
     "e": 10828,
     "st": "+",
     "len": 305,
     "match": 915,
     "num": 915,
     "complete": true
    },
    {
     "prot": "Rv0008c",
     "name": "Rv0008c",
     "s": 11874,
     "e": 12311,
     "st": "-",
     "len": 146,
     "match": 435,
     "num": 438,
     "complete": true
    }
   ],
   "bad": [
    {
     "prot": "G49284_L4_00055-NA",
     "s": 11421,
     "e": 11528,
     "st": "-"
    },
    {
     "prot": "G49284_L4_00060-NA",
     "s": 11555,
     "e": 11692,
     "st": "-"
    }
   ],
   "refMatch": [
    {
     "tag": "L4_CDC1551_00009",
     "match": 216,
     "num": 216
    },
    {
     "tag": "L4_CDC1551_00010",
     "match": 105,
     "num": 105
    },
    {
     "tag": "L4_CDC1551_00011",
     "match": 132,
     "num": 135
    }
   ],
   "amap": [
    {
     "id": "L4_CDC1551_Rv0001",
     "rvs": "Rv0001",
     "s": 1,
     "e": 1524,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_00002",
     "rvs": "Rv0002",
     "s": 2052,
     "e": 3260,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_00003",
     "rvs": "Rv0003",
     "s": 3280,
     "e": 4437,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_Rv0004",
     "rvs": "Rv0004",
     "s": 4434,
     "e": 4997,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_Rv0005",
     "rvs": "Rv0005",
     "s": 5240,
     "e": 7267,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_00006",
     "rvs": "Rv0006",
     "s": 7302,
     "e": 9818,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_Rv0007",
     "rvs": "Rv0007",
     "s": 9914,
     "e": 10828,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_00012",
     "rvs": "Rv0008c",
     "s": 11874,
     "e": 12311,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_Rv0001",
     "rvs": "Rv0001",
     "s": 1,
     "e": 1524,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_00002",
     "rvs": "Rv0002",
     "s": 2052,
     "e": 3260,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_00003",
     "rvs": "Rv0003",
     "s": 3280,
     "e": 4437,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_Rv0004",
     "rvs": "Rv0004",
     "s": 4434,
     "e": 4997,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_Rv0005",
     "rvs": "Rv0005",
     "s": 5240,
     "e": 7267,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_00006",
     "rvs": "Rv0006",
     "s": 7302,
     "e": 9818,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_Rv0007",
     "rvs": "Rv0007",
     "s": 9914,
     "e": 10828,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_00012",
     "rvs": "Rv0008c",
     "s": 11874,
     "e": 12311,
     "st": "-"
    }
   ]
  },
  "case2": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 1980,
   "to": 3275,
   "bakta": [
    {
     "tag": "L4_CDC1551_00002",
     "s": 2052,
     "e": 3260,
     "st": "+",
     "name": "dnaN",
     "product": "DNA polymerase III (beta chain) DnaN (DNA nucleotidyltransferase)",
     "case": 2
    }
   ],
   "ref": [
    {
     "prot": "Rv0002",
     "name": "dnaN",
     "s": 2052,
     "e": 3260,
     "st": "+",
     "len": 403,
     "match": 1209,
     "num": 1209,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551_00002",
     "rvs": "Rv0002",
     "s": 2052,
     "e": 3260,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_00002",
     "rvs": "Rv0002",
     "s": 2052,
     "e": 3260,
     "st": "+"
    }
   ]
  },
  "case3": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 23780,
   "to": 25650,
   "bakta": [
    {
     "tag": "L4_CDC1551_00025",
     "s": 23861,
     "e": 25426,
     "st": "-",
     "name": "fhaA",
     "product": "Conserved protein with FHA domain, FhaA",
     "case": 3
    }
   ],
   "ref": [
    {
     "prot": "Rv0020c",
     "name": "fhaA",
     "s": 23861,
     "e": 25426,
     "st": "-",
     "len": 528,
     "match": 1566,
     "num": 1584,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551_00025",
     "rvs": "Rv0020c",
     "s": 23861,
     "e": 25426,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_00025",
     "rvs": "Rv0020c",
     "s": 23861,
     "e": 25426,
     "st": "-"
    }
   ]
  },
  "case4-start": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 1,
   "to": 1700,
   "bakta": [
    {
     "tag": "L4_CDC1551_00001",
     "s": 34,
     "e": 1524,
     "st": "+",
     "name": "dnaA",
     "product": "Chromosomal replication initiator protein DnaA",
     "case": 4
    }
   ],
   "ref": [
    {
     "prot": "Rv0001",
     "name": "dnaA",
     "s": 1,
     "e": 1524,
     "st": "+",
     "len": 508,
     "match": 1524,
     "num": 1524,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551_Rv0001",
     "rvs": "Rv0001",
     "s": 1,
     "e": 1524,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_Rv0001",
     "rvs": "Rv0001",
     "s": 1,
     "e": 1524,
     "st": "+"
    }
   ]
  },
  "case4-split": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 9850,
   "to": 11000,
   "bakta": [
    {
     "tag": "L4_CDC1551_00007",
     "s": 9963,
     "e": 10160,
     "st": "+",
     "name": "",
     "product": "hypothetical protein",
     "case": 4
    },
    {
     "tag": "L4_CDC1551_00008",
     "s": 10523,
     "e": 10828,
     "st": "+",
     "name": "",
     "product": "hypothetical protein",
     "case": 4
    }
   ],
   "ref": [
    {
     "prot": "Rv0007",
     "name": "Rv0007",
     "s": 9914,
     "e": 10828,
     "st": "+",
     "len": 305,
     "match": 915,
     "num": 915,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551_Rv0007",
     "rvs": "Rv0007",
     "s": 9914,
     "e": 10828,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_Rv0007",
     "rvs": "Rv0007",
     "s": 9914,
     "e": 10828,
     "st": "+"
    }
   ]
  },
  "case4-longer": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 2267650,
   "to": 2272700,
   "bakta": [
    {
     "tag": "L4_CDC1551_02139",
     "s": 2267753,
     "e": 2272573,
     "st": "-",
     "name": "",
     "product": "site-specific DNA-methyltransferase (adenine-specific)",
     "case": 4
    }
   ],
   "ref": [
    {
     "prot": "Rv2024c",
     "name": "Rv2024c",
     "s": 2271056,
     "e": 2272573,
     "st": "-",
     "len": 516,
     "match": 1512,
     "num": 1518,
     "complete": false
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551_02139",
     "rvs": "Rv2024c",
     "s": 2267753,
     "e": 2272573,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_02139",
     "rvs": "Rv2024c",
     "s": 2267753,
     "e": 2272573,
     "st": "-"
    }
   ]
  },
  "case1-kept": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 1979450,
   "to": 1985565,
   "bakta": [
    {
     "tag": "L4_CDC1551_01841",
     "s": 1979541,
     "e": 1980686,
     "st": "-",
     "name": "",
     "product": "GDP-mannose-dependent alpha-(1-6)-phosphatidylinositol dimannoside mannosyltransferase",
     "case": 1
    },
    {
     "tag": "L4_CDC1551_01842",
     "s": 1980912,
     "e": 1982015,
     "st": "+",
     "name": "",
     "product": "Sulfite oxidase-like oxidoreductase",
     "case": 1
    },
    {
     "tag": "L4_CDC1551_01843",
     "s": 1982216,
     "e": 1985053,
     "st": "+",
     "name": "",
     "product": "Membrane transport protein MMPL domain-containing protein",
     "case": 1
    },
    {
     "tag": "L4_CDC1551_01844",
     "s": 1985374,
     "e": 1985562,
     "st": "+",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    }
   ],
   "ref": [],
   "bad": [],
   "refMatch": [
    {
     "tag": "L4_CDC1551_01841",
     "match": 0,
     "num": 0
    },
    {
     "tag": "L4_CDC1551_01842",
     "match": 0,
     "num": 0
    },
    {
     "tag": "L4_CDC1551_01843",
     "match": 0,
     "num": 0
    },
    {
     "tag": "L4_CDC1551_01844",
     "match": 90,
     "num": 150
    }
   ],
   "amap": [
    {
     "id": "L4_CDC1551_01841",
     "rvs": "",
     "s": 1979541,
     "e": 1980686,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_01842",
     "rvs": "",
     "s": 1980912,
     "e": 1982015,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_01843",
     "rvs": "",
     "s": 1982216,
     "e": 1985053,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_01844",
     "rvs": "",
     "s": 1985374,
     "e": 1985562,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_01841",
     "rvs": "",
     "s": 1979541,
     "e": 1980686,
     "st": "-"
    },
    {
     "id": "L4_CDC1551_01842",
     "rvs": "",
     "s": 1980912,
     "e": 1982015,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_01843",
     "rvs": "",
     "s": 1982216,
     "e": 1985053,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_01844",
     "rvs": "",
     "s": 1985374,
     "e": 1985562,
     "st": "+"
    }
   ]
  },
  "case1-dropped": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 11400,
   "to": 11545,
   "bakta": [
    {
     "tag": "L4_CDC1551_00010",
     "s": 11421,
     "e": 11528,
     "st": "-",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    }
   ],
   "ref": [],
   "bad": [
    {
     "prot": "G49284_L4_00055-NA",
     "s": 11421,
     "e": 11528,
     "st": "-"
    }
   ],
   "refMatch": [
    {
     "tag": "L4_CDC1551_00010",
     "match": 105,
     "num": 105
    }
   ],
   "amap": [],
   "corrected": []
  },
  "missing": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 33100,
   "to": 34420,
   "bakta": [
    {
     "tag": "L4_CDC1551_00037",
     "s": 33168,
     "e": 33497,
     "st": "+",
     "name": "",
     "product": "Conserved hypothetical protein",
     "case": 2
    },
    {
     "tag": "L4_CDC1551_00038",
     "s": 34239,
     "e": 36554,
     "st": "+",
     "name": "bioF2",
     "product": "putative 8-amino-7-oxononanoate synthase BioF2 (AONS) (8-amino-7-ketopelargonate synthase) (7-keto-8-amino-pelargonic acid synthetase) (7-KAP synthetase) (L-alanine--pimelyl CoA ligase)",
     "case": 2
    }
   ],
   "ref": [
    {
     "prot": "Rv0030",
     "name": "Rv0030",
     "s": 33168,
     "e": 33497,
     "st": "+",
     "len": 110,
     "match": 330,
     "num": 330,
     "complete": true
    },
    {
     "prot": "Rv0031",
     "name": "Rv0031",
     "s": 33526,
     "e": 33738,
     "st": "+",
     "len": 71,
     "match": 213,
     "num": 213,
     "complete": true
    },
    {
     "prot": "Rv0032",
     "name": "bioF2",
     "s": 34239,
     "e": 36554,
     "st": "+",
     "len": 772,
     "match": 2316,
     "num": 2316,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551_00037",
     "rvs": "Rv0030",
     "s": 33168,
     "e": 33497,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_Rv0031",
     "rvs": "Rv0031",
     "s": 33526,
     "e": 33738,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_00038",
     "rvs": "Rv0032",
     "s": 34239,
     "e": 36554,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_00037",
     "rvs": "Rv0030",
     "s": 33168,
     "e": 33497,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_Rv0031",
     "rvs": "Rv0031",
     "s": 33526,
     "e": 33738,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_00038",
     "rvs": "Rv0032",
     "s": 34239,
     "e": 36554,
     "st": "+"
    }
   ]
  },
  "paralogs": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "kind": "synteny",
   "placements": [
    {
     "prot": "Rv0331",
     "name": "Rv0331",
     "s": 396261,
     "e": 397427,
     "st": "+",
     "refStart": 396201,
     "refEnd": 397367,
     "match": 1164,
     "num": 1167,
     "len": 389,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0332",
     "name": "Rv0332",
     "s": 397502,
     "e": 398287,
     "st": "+",
     "refStart": 397442,
     "refEnd": 398227,
     "match": 786,
     "num": 786,
     "len": 262,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0333",
     "name": "Rv0333",
     "s": 398314,
     "e": 398688,
     "st": "+",
     "refStart": 398254,
     "refEnd": 398628,
     "match": 375,
     "num": 375,
     "len": 125,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0334",
     "name": "rmlA",
     "s": 398718,
     "e": 399584,
     "st": "+",
     "refStart": 398658,
     "refEnd": 399524,
     "match": 867,
     "num": 867,
     "len": 289,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0335c",
     "name": "PE6",
     "s": 399595,
     "e": 400110,
     "st": "-",
     "refStart": 399535,
     "refEnd": 400050,
     "match": 516,
     "num": 516,
     "len": 172,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0336",
     "name": "Rv0336",
     "s": 400252,
     "e": 401763,
     "st": "+",
     "refStart": 400192,
     "refEnd": 401703,
     "match": 1512,
     "num": 1512,
     "len": 504,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0515",
     "name": "Rv0515",
     "s": 400252,
     "e": 401763,
     "st": "+",
     "refStart": 606551,
     "refEnd": 608062,
     "match": 1509,
     "num": 1512,
     "len": 504,
     "complete": true,
     "kept": false
    },
    {
     "prot": "Rv0337c",
     "name": "aspC",
     "s": 401933,
     "e": 403222,
     "st": "-",
     "refStart": 401873,
     "refEnd": 403162,
     "match": 1290,
     "num": 1290,
     "len": 430,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0338c",
     "name": "Rv0338c",
     "s": 403253,
     "e": 405901,
     "st": "-",
     "refStart": 403193,
     "refEnd": 405841,
     "match": 2637,
     "num": 2649,
     "len": 883,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0339c",
     "name": "Rv0339c",
     "s": 406010,
     "e": 408508,
     "st": "-",
     "refStart": 405950,
     "refEnd": 408448,
     "match": 2499,
     "num": 2499,
     "len": 833,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0340",
     "name": "Rv0340",
     "s": 408694,
     "e": 409233,
     "st": "+",
     "refStart": 408634,
     "refEnd": 409173,
     "match": 540,
     "num": 540,
     "len": 180,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0341",
     "name": "iniB",
     "s": 409422,
     "e": 410861,
     "st": "+",
     "refStart": 409362,
     "refEnd": 410801,
     "match": 1431,
     "num": 1440,
     "len": 480,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0510",
     "name": "hemC",
     "s": 603224,
     "e": 604153,
     "st": "+",
     "refStart": 601857,
     "refEnd": 602786,
     "match": 930,
     "num": 930,
     "len": 310,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0511",
     "name": "hemD",
     "s": 604186,
     "e": 605883,
     "st": "+",
     "refStart": 602819,
     "refEnd": 604516,
     "match": 1698,
     "num": 1698,
     "len": 566,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0512",
     "name": "hemB",
     "s": 605969,
     "e": 606958,
     "st": "+",
     "refStart": 604602,
     "refEnd": 605591,
     "match": 990,
     "num": 990,
     "len": 330,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0513",
     "name": "Rv0513",
     "s": 606971,
     "e": 607519,
     "st": "+",
     "refStart": 605604,
     "refEnd": 606152,
     "match": 549,
     "num": 549,
     "len": 183,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0514",
     "name": "Rv0514",
     "s": 607516,
     "e": 607815,
     "st": "+",
     "refStart": 606149,
     "refEnd": 606448,
     "match": 300,
     "num": 300,
     "len": 100,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0336",
     "name": "Rv0336",
     "s": 607918,
     "e": 609429,
     "st": "+",
     "refStart": 400192,
     "refEnd": 401703,
     "match": 1509,
     "num": 1512,
     "len": 504,
     "complete": true,
     "kept": false
    },
    {
     "prot": "Rv0515",
     "name": "Rv0515",
     "s": 607918,
     "e": 609429,
     "st": "+",
     "refStart": 606551,
     "refEnd": 608062,
     "match": 1512,
     "num": 1512,
     "len": 504,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0516c",
     "name": "Rv0516c",
     "s": 609426,
     "e": 609902,
     "st": "-",
     "refStart": 608059,
     "refEnd": 608535,
     "match": 477,
     "num": 477,
     "len": 159,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0517",
     "name": "Rv0517",
     "s": 610113,
     "e": 611423,
     "st": "+",
     "refStart": 608746,
     "refEnd": 610056,
     "match": 1311,
     "num": 1311,
     "len": 437,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0518",
     "name": "Rv0518",
     "s": 611555,
     "e": 612250,
     "st": "+",
     "refStart": 610188,
     "refEnd": 610883,
     "match": 696,
     "num": 696,
     "len": 232,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0519c",
     "name": "Rv0519c",
     "s": 612539,
     "e": 613441,
     "st": "-",
     "refStart": 611172,
     "refEnd": 612074,
     "match": 903,
     "num": 903,
     "len": 301,
     "complete": true,
     "kept": true
    },
    {
     "prot": "Rv0520",
     "name": "Rv0520",
     "s": 613622,
     "e": 613972,
     "st": "+",
     "refStart": 612255,
     "refEnd": 612605,
     "match": 351,
     "num": 351,
     "len": 117,
     "complete": true,
     "kept": true
    }
   ]
  },
  "paralogs-1": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 400170,
   "to": 402300,
   "bakta": [
    {
     "tag": "L4_CDC1551_00357",
     "s": 400321,
     "e": 401763,
     "st": "+",
     "name": "",
     "product": "Conserved 13E12 repeat family protein",
     "case": 4
    },
    {
     "tag": "L4_CDC1551_00358",
     "s": 401933,
     "e": 403222,
     "st": "-",
     "name": "aspC",
     "product": "putative aspartate aminotransferase AspC (transaminase A) (ASPAT)",
     "case": 2
    }
   ],
   "ref": [
    {
     "prot": "Rv0336",
     "name": "Rv0336",
     "s": 400252,
     "e": 401763,
     "st": "+",
     "len": 504,
     "match": 1512,
     "num": 1512,
     "complete": true
    },
    {
     "prot": "Rv0337c",
     "name": "aspC",
     "s": 401933,
     "e": 403222,
     "st": "-",
     "len": 430,
     "match": 1290,
     "num": 1290,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551_Rv0336",
     "rvs": "Rv0336",
     "s": 400252,
     "e": 401763,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_00358",
     "rvs": "Rv0337c",
     "s": 401933,
     "e": 403222,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "placed": [
    {
     "prot": "Rv0336",
     "s": 400252,
     "e": 401763,
     "st": "+",
     "match": 1512,
     "num": 1512
    },
    {
     "prot": "Rv0337c",
     "s": 401933,
     "e": 403222,
     "st": "-",
     "match": 1290,
     "num": 1290
    },
    {
     "prot": "Rv0515",
     "s": 400252,
     "e": 401763,
     "st": "+",
     "match": 1509,
     "num": 1512
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_Rv0336",
     "rvs": "Rv0336",
     "s": 400252,
     "e": 401763,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_00358",
     "rvs": "Rv0337c",
     "s": 401933,
     "e": 403222,
     "st": "-"
    }
   ]
  },
  "paralogs-2": {
   "genome": "M. tuberculosis CDC1551",
   "base": "L4_CDC1551",
   "contig": "L4_CDC1551",
   "from": 607300,
   "to": 609900,
   "bakta": [
    {
     "tag": "L4_CDC1551_00543",
     "s": 606971,
     "e": 607519,
     "st": "+",
     "name": "",
     "product": "putative conserved transmembrane protein",
     "case": 2
    },
    {
     "tag": "L4_CDC1551_00544",
     "s": 607516,
     "e": 607815,
     "st": "+",
     "name": "",
     "product": "putative transmembrane protein",
     "case": 2
    },
    {
     "tag": "L4_CDC1551_00545",
     "s": 607987,
     "e": 609429,
     "st": "+",
     "name": "",
     "product": "Conserved 13E12 repeat family protein",
     "case": 4
    },
    {
     "tag": "L4_CDC1551_00546",
     "s": 609426,
     "e": 609902,
     "st": "-",
     "name": "",
     "product": "putative anti-anti-sigma factor",
     "case": 2
    }
   ],
   "ref": [
    {
     "prot": "Rv0513",
     "name": "Rv0513",
     "s": 606971,
     "e": 607519,
     "st": "+",
     "len": 183,
     "match": 549,
     "num": 549,
     "complete": true
    },
    {
     "prot": "Rv0514",
     "name": "Rv0514",
     "s": 607516,
     "e": 607815,
     "st": "+",
     "len": 100,
     "match": 300,
     "num": 300,
     "complete": true
    },
    {
     "prot": "Rv0515",
     "name": "Rv0515",
     "s": 607918,
     "e": 609429,
     "st": "+",
     "len": 504,
     "match": 1512,
     "num": 1512,
     "complete": true
    },
    {
     "prot": "Rv0516c",
     "name": "Rv0516c",
     "s": 609426,
     "e": 609902,
     "st": "-",
     "len": 159,
     "match": 477,
     "num": 477,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551_00543",
     "rvs": "Rv0513",
     "s": 606971,
     "e": 607519,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_00544",
     "rvs": "Rv0514",
     "s": 607516,
     "e": 607815,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_Rv0515",
     "rvs": "Rv0515",
     "s": 607918,
     "e": 609429,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551_00546",
     "rvs": "Rv0516c",
     "s": 609426,
     "e": 609902,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "placed": [
    {
     "prot": "Rv0336",
     "s": 607918,
     "e": 609429,
     "st": "+",
     "match": 1509,
     "num": 1512
    },
    {
     "prot": "Rv0513",
     "s": 606971,
     "e": 607519,
     "st": "+",
     "match": 549,
     "num": 549
    },
    {
     "prot": "Rv0514",
     "s": 607516,
     "e": 607815,
     "st": "+",
     "match": 300,
     "num": 300
    },
    {
     "prot": "Rv0515",
     "s": 607918,
     "e": 609429,
     "st": "+",
     "match": 1512,
     "num": 1512
    },
    {
     "prot": "Rv0516c",
     "s": 609426,
     "e": 609902,
     "st": "-",
     "match": 477,
     "num": 477
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551_00543",
     "rvs": "Rv0513",
     "s": 606971,
     "e": 607519,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_00544",
     "rvs": "Rv0514",
     "s": 607516,
     "e": 607815,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_Rv0515",
     "rvs": "Rv0515",
     "s": 607918,
     "e": 609429,
     "st": "+"
    },
    {
     "id": "L4_CDC1551_00546",
     "rvs": "Rv0516c",
     "s": 609426,
     "e": 609902,
     "st": "-"
    }
   ]
  },
  "small-protein": {
   "genome": "E. coli K-12",
   "base": "K12self",
   "contig": "K12self",
   "from": 3056150,
   "to": 3058450,
   "bakta": [
    {
     "tag": "K12self_03011",
     "s": 3056187,
     "e": 3056789,
     "st": "+",
     "name": "fau",
     "product": "5-formyltetrahydrofolate cyclo-ligase",
     "case": 4
    },
    {
     "tag": "K12self_03013",
     "s": 3057178,
     "e": 3058410,
     "st": "-",
     "name": "serA",
     "product": "phosphoglycerate dehydrogenase",
     "case": 2
    }
   ],
   "ref": [
    {
     "prot": "b2912",
     "name": "fau",
     "s": 3056241,
     "e": 3056789,
     "st": "+",
     "len": 183,
     "match": 549,
     "num": 549,
     "complete": true
    },
    {
     "prot": "b2913",
     "name": "serA",
     "s": 3057178,
     "e": 3058410,
     "st": "-",
     "len": 411,
     "match": 1233,
     "num": 1233,
     "complete": true
    },
    {
     "prot": "b4665",
     "name": "ibsC",
     "s": 3056890,
     "e": 3056949,
     "st": "-",
     "len": 20,
     "match": 60,
     "num": 60,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "K12self_b2912",
     "rvs": "b2912",
     "s": 3056241,
     "e": 3056789,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "K12self_b4665",
     "rvs": "b4665",
     "s": 3056890,
     "e": 3056949,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "K12self_03013",
     "rvs": "b2913",
     "s": 3057178,
     "e": 3058410,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "miniprot": [
    {
     "prot": "b2912",
     "s": 3056241,
     "e": 3056789,
     "st": "+",
     "match": 549,
     "num": 549
    },
    {
     "prot": "b2913",
     "s": 3057178,
     "e": 3058410,
     "st": "-",
     "match": 1233,
     "num": 1233
    }
   ],
   "corrected": [
    {
     "id": "K12self_b2912",
     "rvs": "b2912",
     "s": 3056241,
     "e": 3056789,
     "st": "+"
    },
    {
     "id": "K12self_b4665",
     "rvs": "b4665",
     "s": 3056890,
     "e": 3056949,
     "st": "-"
    },
    {
     "id": "K12self_03013",
     "rvs": "b2913",
     "s": 3057178,
     "e": 3058410,
     "st": "-"
    }
   ],
   "exact": {
    "prot": "b4665",
    "protein": "MMRLVIILIVLLLISFSAY*",
    "offset": 3056859,
    "dna": "GGAGGATTGCTCCTCCCCTGAGACTGACTGTTAATAAGCGCTGAAACTTATGAGTAACAGTACAATCAGTATGATGACAAGTCGCATCATAACCCTTCTCCTTCAAGCCCTCGCTTCGGT",
    "elsewhere": [
     {
      "s": 3194766,
      "e": 3194825,
      "st": "+",
      "match": 51,
      "num": 60
     }
    ]
   }
  },
  "contig-end": {
   "genome": "M. tuberculosis CDC1551 in 40 contigs",
   "base": "L4_CDC1551d",
   "contig": "contig_10",
   "from": 100300,
   "to": 102146,
   "bakta": [
    {
     "tag": "L4_CDC1551d_01398",
     "s": 100318,
     "e": 101307,
     "st": "-",
     "name": "mdh",
     "product": "putative malate dehydrogenase Mdh",
     "case": 2
    },
    {
     "tag": "L4_CDC1551d_01399",
     "s": 101451,
     "e": 101756,
     "st": "-",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    },
    {
     "tag": "L4_CDC1551d_01400",
     "s": 101760,
     "e": 101957,
     "st": "-",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    }
   ],
   "ref": [
    {
     "prot": "Rv1239c",
     "name": "corA",
     "s": 101478,
     "e": 102146,
     "st": "+",
     "len": 367,
     "match": 666,
     "num": 669,
     "complete": false
    },
    {
     "prot": "Rv1240",
     "name": "mdh",
     "s": 100318,
     "e": 101307,
     "st": "-",
     "len": 330,
     "match": 990,
     "num": 990,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [
    {
     "tag": "L4_CDC1551d_01399",
     "match": 300,
     "num": 303
    },
    {
     "tag": "L4_CDC1551d_01400",
     "match": 192,
     "num": 195
    }
   ],
   "amap": [
    {
     "id": "L4_CDC1551d_01398",
     "rvs": "Rv1240",
     "s": 100318,
     "e": 101307,
     "st": "-",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551d_Rv1239c_2",
     "rvs": "Rv1239c",
     "s": 101478,
     "e": 102146,
     "st": "+",
     "pseudo": false,
     "partial": "3'",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551d_01398",
     "rvs": "Rv1240",
     "s": 100318,
     "e": 101307,
     "st": "-"
    },
    {
     "id": "L4_CDC1551d_Rv1239c",
     "rvs": "Rv1239c",
     "s": 101478,
     "e": 102146,
     "st": "+"
    }
   ],
   "translate": {
    "offset": 99099,
    "dna": "GGCGGCTTTGGGGGCAGCGGGGGTCGCGGCGGGCTACTGCTCGGCACGGGCGGCGCCGGCGGCGACGGCGGCGTCGGGGGCGACGGGGGCATCGGCGCTCAAGGCGGCAGCGGCGGCAACGGCGGCAACGGCGGGATCGGCGCCGACGGCATGGCCAACCAGGACGGCGACGGCGGTGACGGCGGCAACGGCGGCGACGGCGGGGCCGGCGGGGCCGGTGGCGTCGGCGGAAACGGCGGGACCGGCGGTGCGGGCGGACTGTTCGGACAGTCGGGCAGCCCCGGCTCCGGCGCGGCCGGGGGCCTCGGCGGCGCGGGCGGCAACGGCGGCGCGGGCGGCGGCGGCGGGACCGGGTTCAACCCCGGCGCCCCCGGCGATCCCGGTACTCAAGGCGCTACCGGCGCCAACGGTCAGCACGGCCTGAACGGCTGACCCGAGGGCAGAGCGGCAGTCATTATTCGCGCAACGGGTCGGTCCACTTCAGATCGGCGAATCGGGCAAAGTCGGTGTCACTGGAATGGATCTCGGCGTCGTATTCGATGCCGTAGGCGGCCAGTTGCACATCGGTGGTTAGGTGGCTGGCTGTGCCGAGCTTGTCGAGCAGGCCCAACGCGATGTCCAGGTGGCGCGGACCCGCCGTGAGTAGGTCCACGTTCGGCTGCGAAAGCCACTCGCGCACATAGGCGATCGCATCCGCGGTTGGCAGTGGCGCGGCGAGCACGCGGGCGCTGGTGGCGATCCGTAGGAACCCGAACAACGCCGGATACGTCAGCCCGATACGGGTGTGGCCGTTGACGGTGTCTTGCCACCACGCATGCGCGCGCCGGTGCTGCGGGAATCCGGTGATGACCGCGTAGAGCAGCAGATTGATGTCAGGGATGATCACCGGGCCCGACGCGTGGCATCCAGCAGCGCCTCATCCTCCAGTTCGTCGGCCAACTTGTTGAAGCCGGCCAGATCCAACCCGGAACGCACAGCCGACTCATGCGGCTCCAACCGATACTGCTCCTGCCGTTTCACCGGCGGCGCCAGCGCTCTGCGCAGCGCATCGTTGATGACCTGCTTCATCGGGCGGCGTTCGCGATGCACTGCGTCTTCGACCAGCCGGACGACGTCGTCATCGAGCGTCAAGGTGGTGCGCATACAGGTCATCATGATGCCTGGACAAACACACGTCAAGATGACTCGGGCTCCGCTCAGTGCGCGGCTGACCTGCGCTCAGATCAGGCCGAGCTCGGTGACCGCGCTGCGCTCGTCAGCCAACTCGGCGGTTGACTTGTCGATCCGGCCGCGGGAGAACTCGTCGATCTCCAAGCCGCTCACGATCGTCCAGTTGCCGCCCTTGGTGGTGACCGGAAACGAGGAGATCAAGCCCTCCGGCACCCCGTAGGACCCGTCGGAGACGACGGCCATCGAGACCCAATCGTCCGCCGGCGTCCCCAGCAACCAGTCCCGGGCAGCGTCGATGGTTGCCGACGCGGCCGAGGCGGCCGACGACGCGCCGCGCGCATCGATGATCGCCGCACCGCGCTTGGCGACCGTCGGGATGAATTCATCCTCGATCCAGGCCTGGTCGTTGACCACTTCGGCCGCGTTCTTTCCGGCGACCTCCGCGTGGAACAGGTCGGGGTACTGGGTGGCCGAGTGATTGCCCCAGATCGTCATCTTCTTGATGTCGGTGACCGCCGCGCCGGTCTTGGCGGCCAGCTGCGAGATCGCCCGATTGTGGTCCAGCCGGGTGAGCGCCGAGAACCGCTCGCGGGGAATGTCGGGCGCATTGGTCATCGCGATCAGCGCGTTGGTGTTGGCGGGGTTGCCGGTCACCCCGACGCGAACGTCATCCGCGGCGACAGCGTTGAGGGCTTTGCCCTGAGCGGTGAAGATCGCGCCGTTGGCCTCCAGCAGGTCACTTCGCTCCATGCCCGCGCCCCGGGGGCGGGCTCCGACCAGCAGGGCCAGGCTCACACCATCGAAGATCTTCTGGGGATCGGAACCGATCTCCACCCCGGACAACAGCGGGAAAGCGCAGTCGTCGAGTTCCATCACCACACCCTCGAGCGCCTGCAGTGCCGGCTCGATCTCGAGCAGCCGCAGCTCGATCGGACGGTCAGGGCCCAGCAAAGAGCCGCTGGCCAGGCGGAACAACAGGCTGTAGCCGATTTGGCCGGCGGCGCCGGTAACGGCGACCTTGAGAGGACTAGCGCTCACGTCGGTTCTGCTCCTTACGGACTTGGCTTCGGGTTGAAACTAGCGCACCGGCACCGGACGGACATCTGGGGTCCAACCCCCACAGCTACGACTGACGTTTACCAAGTTCCATACCGGGACGCGTAACCTGGAGCACCGCCCGATCAGATCCACCGCGATCGGAGGTTGACGTGTTCCCAGGGTTTGACGCATTGCCCGAAGTGCTGCGACCGGTCGCGCGACCCCAGCCGCCGAACGCACACCCCGTTGCCCAGCCACCGGCCCAAGCCTTGGTCGACTGCGGTGTCTACGTCTGCGGCCAGCGACTGCCCGGCAAGTACACCTACGCCGCCGCGCTGCGCGAGGTGCGCGAGATCGAACTGACCGGGCAGGAGGCGTTCGTCTGGATCGGGCTGCACGAGCCCGATGAAAACCAGATGCAGGACGTAGCAGACGTTTTCGGGTTGCACCCGTTAGCCGTTGAGGACGCCGTGCACGCGCACCAGCGACCCAAGTTGGAGCGCTACGACGAGACGCTGTTCCTCGTCCTCAAGACCGTCAACTACGTCCCGCACGAATCGGTGGTACTGGCCCGCGAGATCGTCGAAACCGGCGAGATCATGATCTTCGTCGGCAAGGATTTCGTGGTCACCGTCCGCCACGGCGAACACGGCGGGTTATCCGAGGTGCGTAAGCGGATGGATGCCGACCCCGAACATTTGCGGTTGGGACCGTATGCGGTGATGCACGCGATCGCCGACTACGTGGTCGACCACTACCTCGAGGTGACCAATCTCATGGAGACCGATATCGACAGCATCGAGGAAGTAGCGTTCGCGCCGGGCCGCAAGCTCGACATC",
    "contigLength": 102146,
    "proteins": {
     "Rv1240": [
      "VSASPLKVAVTGAAGQIGYSLLFRLASGSLLGPDRPIELRLLEIEPALQALEGVVMELDDCAFPLLSGVEIGSDPQKIFDGVSLALLVGARPRGAGMERSDLLEANGAIFTAQGKALNAVAADDVRVGVTGNPANTNALIAMTNAPDIPRERFSALTRLDHNRAISQLAAKTGAAVTDIKKMTIWGNHSATQYPDLFHAEVAGKNAAEVVNDQAWIEDEFIPTVAKRGAAIIDARGASSAASAASATIDAARDWLLGTPADDWVSMAVVSDGSYGVPEGLISSFPVTTKGGNWTIVSGLEIDEFSRGRIDKSTAELADERSAVTELGLI*",
      990
     ],
     "Rv1239c": [
      "VFPGFDALPEVLRPVARPQPPNAHPVAQPPAQALVDCGVYVCGQRLPGKYTYAAALREVREIELTGQEAFVWIGLHEPDENQMQDVADVFGLHPLAVEDAVHAHQRPKLERYDETLFLVLKTVNYVPHESVVLAREIVKTGEIMIFVGKDFVVTVRHGEHGGLSEVRKRMDADPEHLRLGPYAVMHAIADYVVDHYLEVTNLMETDIDSIEEVAFAPGRKLDIEPIYLLKREVVELRRCVNPLSTAFQRMQTESKDLISKEVRRYLRDVADHQTEAADQIASYDDMLNSLVQAALARVGMQQNMDMRKISAWAGIIAVPTMIAGIYGMNFHFMPELDSRWGYPTVIGGMVLICLFLYHVFRNRNWL*",
      1101
     ]
    }
   }
  },
  "contig-end-2": {
   "genome": "M. tuberculosis CDC1551 in 40 contigs",
   "base": "L4_CDC1551d",
   "contig": "contig_9",
   "from": 3700,
   "to": 4591,
   "bakta": [
    {
     "tag": "L4_CDC1551d_01310",
     "s": 2902,
     "e": 4083,
     "st": "+",
     "name": "sugC",
     "product": "putative sugar-transport ATP-binding protein ABC transporter SugC",
     "case": 2
    },
    {
     "tag": "L4_CDC1551d_01311",
     "s": 4160,
     "e": 4513,
     "st": "-",
     "name": "",
     "product": "hypothetical protein",
     "case": 4
    }
   ],
   "ref": [
    {
     "prot": "Rv1238",
     "name": "sugC",
     "s": 2902,
     "e": 4083,
     "st": "+",
     "len": 394,
     "match": 1182,
     "num": 1182,
     "complete": true
    },
    {
     "prot": "Rv1239c",
     "name": "corA",
     "s": 4160,
     "e": 4591,
     "st": "-",
     "len": 367,
     "match": 432,
     "num": 432,
     "complete": false,
     "piece": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "L4_CDC1551d_01310",
     "rvs": "Rv1238",
     "s": 2902,
     "e": 4083,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "L4_CDC1551d_Rv1239c",
     "rvs": "Rv1239c",
     "s": 4160,
     "e": 4591,
     "st": "-",
     "pseudo": false,
     "partial": "5'",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "L4_CDC1551d_01310",
     "rvs": "Rv1238",
     "s": 2902,
     "e": 4083,
     "st": "+"
    },
    {
     "id": "L4_CDC1551d_Rv1239c",
     "rvs": "Rv1239c",
     "s": 4160,
     "e": 4591,
     "st": "-"
    }
   ],
   "translate": {
    "offset": 1701,
    "dna": "GCTTTTGGCCGGGTTGGCGCTGGTCCCCGAGGATCTGCTAAGAGCAGCGCAGGTTGACGGCGCCAGCGCCTGGCGGCGGTTGACGAAGGTCATCTTGCCGATGATCAAGCCGGCGATCGTGGTTGCTCTGCTCTTCAGGACCCTGGACGCTTTCCGGATTTTCGACAACATCTATGTGCTGACCGGCGGCAGCAACAACACCGGATCGGTGTCGATCTTGGGCTACGACAACCTGTTCAAGGGGTTCAACGTGGGCCTTGGTTCGGCGATCAGCGTGCTGATCTTTGGCTGCGTGGCCGTCATTGCGTTCATTTTCATCAAGTTGTTCGGCGCCGCGGCGCCCGGGGGTGAGCCAAGTGGGCGTTGAACGGGTGGGCGCGCGGCGCGCCACGTATTGGGCCGTCCTGGACACTTTGGTCGTGGGGTACGCGTTGCTCCCGGTGCTGTGGATTTTCAGCCTGTCACTCAAGCCGACGTCAACGGTCAAGGACGGCAAGCTGATTCCGTCGACGGTGACTTTCGACAACTATCGTGGCATCTTCCGGGGCGACTTGTTCAGCTCAGCGCTGATCAACTCCATCGGAATCGGCCTGATCACCACCGTGATCGCGGTGGTGCTCGGCGCGATGGCGGCCTACGCGGTTGCCCGGCTGGAATTTCCGGGCAAGCGGCTGCTAATCGGGGCTGCCTTGCTGATCACGATGTTCCCGTCGATCTCTTTGGTCACACCATTGTTCAACATCGAACGTGCCATCGGCCTGTTCGACACCTGGCCGGGGTTGATCTTGCCGTACATCACCTTCGCGTTGCCGCTCGCGATCTACACCCTGTCGGCGTTCTTCCGGGAGATCCCTTGGGATCTGGAAAAGGCGGCCAAGATGGACGGTGCAACGCCCGGTCAGGCTTTCCGGAAGGTGATCGTACCGCTGGCGGCGCCGGGCTTGGTGACCGCTGCAATCCTGGTGTTCATTTTCGCCTGGAACGATCTGCTGCTCGCGTTGTCGCTGACCGCTACCAAGGCGGCGATTACCGCGCCGGTGGCCATCGCCAACTTCACCGGCAGTTCGCAATTCGAGGAGCCGACCGGCTCGATCGCGGCCGGCGCGATCGTGATTACGATCCCGATCATCGTCTTTGTTTTAATCTTCCAACGACGGATTGTCGCCGGGTTGACCTCTGGCGCTGTGAAGGGATAGCGCGATGGCCGAGATTGTGTTGGACCACGTCAACAAGAGTTACCCCGACGGTCACACAGCGGTGCGCGACCTCAACCTCACCATCGCCGACGGCGAATTTCTGATCCTGGTAGGGCCTTCCGGTTGTGGCAAGACCACGACGCTGAATATGATTGCTGGGCTTGAAGATATCTCGTCGGGAGAACTGCGCATCGCCGGTGAGCGGGTAAACGAGAAGGCGCCAAAGGACCGTGACATCGCGATGGTGTTCCAGTCGTACGCGCTTTACCCGCATATGACGGTGCGCCAGAACATCGCGTTCCCGCTGACCCTGGCGAAGATGAGAAAGGCCGACATCGCGCAGAAGGTCTCCGAGACTGCAAAAATCCTTGACCTGACCAACCTTCTGGATCGCAAGCCCTCACAATTGTCGGGTGGTCAGCGACAGCGGGTCGCGATGGGCAGGGCAATCGTGCGCCATCCCAAAGCATTCCTGATGGACGAGCCGCTGTCGAACTTGGACGCGAAGTTGCGGGTCCAGATGCGCGGCGAGATTGCCCAGCTGCAGCGGAGGCTGGGTACCACCACCGTCTACGTCACCCACGACCAGACCGAGGCAATGACGCTGGGCGATCGCGTGGTAGTGATGTACGGGGGCATCGCACAGCAGATCGGCACCCCTGAGGAGCTTTACGAACGGCCCGCCAATCTGTTTGTCGCGGGCTTTATCGGCTCGCCGGCCATGAATTTCTTCCCTGCCAGGCTGACCGCGATCGGACTGACCCTGCCGTTCGGTGAGGTGACGCTGGCCCCCGAAGTCCAGGGGGTGATCGCAGCGCACCCGAAACCGGAAAACGTCATCGTAGGCGTGCGGCCGGAGCATATCCAGGACGCAGCATTGATCGACGCGTATCAACGCATCAGGGCGCTGACCTTCCAGGTGAAGGTCAACTTGGTCGAGTCTTTAGGCGCCGACAAATATCTGTATTTCACTACCGAGAGCCCGGCTGTGCACTCGGTTCAGTTGGACGAGTTGGCGGAGGTAGAGGGGGAGTCGGCGTTACACGAAAATCAGTTCGTGGCAAGGGTTCCCGCCGAGTCCAAGGTAGCCATCGGGCAGTCGGTCGAGTTGGCTTTCGATACCGCCAGACTTGCCGTCTTCGACGCCGACTCCGGTGCGAACCTGACCATTCCGCACCGCGCCTAATGGCGGCGAGCGGACACATAAGCCCCCGCCACGCCGAAGGATTTGGAGCTTTTTGCGTCTGTTCGCCGACGCGAAGCTAGAGCCAGTTTCTGTTGCGGAAGACGTGGTAGAGGAACAGACAGATAAGGACCATCCCGCCGATCACTGTCGGGTAACCCCACCTGGAGTCCAGCTCGGGCATGAAGTGAAAGTTCATGCCATAGATGCCCGCGATCATGGTGGGGACCGCGATGATACCTGCCCACGCGGATATCTTGCGCATGTCCATGTTTTGCTGCATGCCGACCCGGGCGAGCGCGGCCTGCACCAGCGAGTTGAGCATGTCGTCGTAGCTGGCGATCTGGTCGGCGGCCTCGGTCTGGTGGTCGGCGACGTCGCGCAGGTAGCGCCGCACTTCTTTCGAAATGAGGTCTTTGCTCTCGGTCTGCATGCGCTGGAATGCGGTCGATAGCGGATTCACGCACCGGCGCAACTCGACCACTTCCCGCTTGAGCAGATAGATCGGTTC",
    "contigLength": 4591,
    "proteins": {
     "Rv1238": [
      "MAEIVLDHVNKSYPDGHTAVRDLNLTIADGEFLILVGPSGCGKTTTLNMIAGLEDISSGELRIAGERVNEKAPKDRDIAMVFQSYALYPHMTVRQNIAFPLTLAKMRKADIAQKVSETAKILDLTNLLDRKPSQLSGGQRQRVAMGRAIVRHPKAFLMDEPLSNLDAKLRVQMRGEIAQLQRRLGTTTVYVTHDQTEAMTLGDRVVVMYGGIAQQIGTPEELYERPANLFVAGFIGSPAMNFFPARLTAIGLTLPFGEVTLAPEVQGVIAAHPKPENVIVGVRPEHIQDAALIDAYQRIRALTFQVKVNLVESLGADKYLYFTTESPAVHSVQLDELAEVEGESALHENQFVARVPAESKVAIGQSVELAFDTARLAVFDADSGANLTIPHRA*",
      1182
     ],
     "Rv1239c": [
      "VFPGFDALPEVLRPVARPQPPNAHPVAQPPAQALVDCGVYVCGQRLPGKYTYAAALREVREIELTGQEAFVWIGLHEPDENQMQDVADVFGLHPLAVEDAVHAHQRPKLERYDETLFLVLKTVNYVPHESVVLAREIVKTGEIMIFVGKDFVVTVRHGEHGGLSEVRKRMDADPEHLRLGPYAVMHAIADYVVDHYLEVTNLMETDIDSIEEVAFAPGRKLDIEPIYLLKREVVELRRCVNPLSTAFQRMQTESKDLISKEVRRYLRDVADHQTEAADQIASYDDMLNSLVQAALARVGMQQNMDMRKISAWAGIIAVPTMIAGIYGMNFHFMPELDSRWGYPTVIGGMVLICLFLYHVFRNRNWL*",
      1101
     ]
    }
   }
  },
  "frameshift": {
   "genome": "M. bovis AF2122/97",
   "base": "La1_Mbovis",
   "contig": "La1_Mbovis",
   "from": 688900,
   "to": 690500,
   "bakta": [
    {
     "tag": "La1_Mbovis_00621",
     "s": 688064,
     "e": 689278,
     "st": "+",
     "name": "mce2A",
     "product": "Mce-family protein Mce2A",
     "case": 2
    },
    {
     "tag": "La1_Mbovis_00622",
     "s": 689275,
     "e": 690306,
     "st": "+",
     "name": "",
     "product": "Mce-family protein Mce1B",
     "case": 4
    },
    {
     "tag": "La1_Mbovis_00623",
     "s": 690303,
     "e": 691748,
     "st": "+",
     "name": "mce2C",
     "product": "Mce-family protein Mce2C",
     "case": 2
    }
   ],
   "ref": [
    {
     "prot": "Rv0589",
     "name": "mce2A",
     "s": 688064,
     "e": 689278,
     "st": "+",
     "len": 405,
     "match": 1209,
     "num": 1215,
     "complete": true
    },
    {
     "prot": "Rv0590",
     "name": "mce2B",
     "s": 689275,
     "e": 690103,
     "st": "+",
     "len": 276,
     "match": 828,
     "num": 829,
     "complete": true
    },
    {
     "prot": "Rv0590A",
     "name": "Rv0590A",
     "s": 690052,
     "e": 690306,
     "st": "+",
     "len": 85,
     "match": 255,
     "num": 255,
     "complete": true
    },
    {
     "prot": "Rv0591",
     "name": "mce2C",
     "s": 690303,
     "e": 691748,
     "st": "+",
     "len": 482,
     "match": 1440,
     "num": 1446,
     "complete": true
    }
   ],
   "bad": [],
   "refMatch": [],
   "amap": [
    {
     "id": "La1_Mbovis_00621",
     "rvs": "Rv0589",
     "s": 688064,
     "e": 689278,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "La1_Mbovis_Rv0590",
     "rvs": "Rv0590",
     "s": 689275,
     "e": 690103,
     "st": "+",
     "pseudo": true,
     "partial": "",
     "note": "frameshift or lost stop codon: its frame runs on 67 residues past the end of the reference protein"
    },
    {
     "id": "La1_Mbovis_Rv0590A",
     "rvs": "Rv0590A",
     "s": 690052,
     "e": 690306,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "La1_Mbovis_00623",
     "rvs": "Rv0591",
     "s": 690303,
     "e": 691748,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "La1_Mbovis_00621",
     "rvs": "Rv0589",
     "s": 688064,
     "e": 689278,
     "st": "+"
    },
    {
     "id": "La1_Mbovis_Rv0590",
     "rvs": "Rv0590",
     "s": 689275,
     "e": 690103,
     "st": "+"
    },
    {
     "id": "La1_Mbovis_Rv0590A",
     "rvs": "Rv0590A",
     "s": 690052,
     "e": 690306,
     "st": "+"
    },
    {
     "id": "La1_Mbovis_00623",
     "rvs": "Rv0591",
     "s": 690303,
     "e": 691748,
     "st": "+"
    }
   ],
   "translate": {
    "offset": 686863,
    "dna": "ATAACCATTGGCCTGGTTGGTGGTTTCGTCTTCAGTGTCTTCATCCAACACGTCTCGGCCGGCGCCTACGTGGGCACGCTCACCTTGGTCACCGGTCTACCCGAGGTGATCATCTCGGTGGTCAAGTCGGCGACGTTCGGCCTGATCGCTGGCCTAGTCGGCTGTTACCGCGGGCTGACCACGAAAGGCGGCCCCAAGGGAGTTGGAACCGCCGTCAACGAAACCCTGGTGCTGTGCGTGATCGCGCTGTTCGCGACCAATGTGGTGTTGACCACGATCGGCGTGCGGTTCGGGACGGGACACTAGCATGGTGGAGTCTTCAACGGCATCAGCGGCAGCCGTATTGCGGGCCCGCTACCCACGCACAGCCGCCAGCCTTGACCGCTACGGCGGCGGCACGGCCCGAAGACTTGAGCGGACAGGGACTTTCGCGAGATTCACCCGGATCAGCGTCGTGCAGATCGGCTGGGCACTGCGTCGCTATCGCCGGGAGACGCTGCGCCTGGTCGCCGAGATCGGGATGGGCACCGGCGCGATGGCCGTCGTCGGCGGCACGGTCGCGATCATCGGTTTTGTGACGCTGTCCGGCGGCTCGCTGATCGCCATCCAGGGCTTCGCGTCGCTGGGCAACATCGGTGTCGAGGCGTTTACCGGATTCTTTGCCGCACTGGCCAACACACGCGTCGCTGCGCCCATTGTCTCCGGTGTCGCGCTGGCCGCGACGGTGGGCGCCGGCGCCACCGCACAGTTAGGTGCCATGCGGATCAGTGAGGAGATCGACGCGCTGGAAGTGATGGGCATCAAGTCGATTTCGTTTCTGGTCTCCACTCGGATTCTAGGAGGGCTGGTGGTGATCATGCCGCTGTACGCGCTCGCTCTCGACATGGCTTTCACCTCTGGTCAGGTGGTCACAACCGTGTTCTACGGCCAGTCCAACGGCACCTATGAGCACTACTTCCGCACCTTCCTGCGCCCAGAGGATGTGGGTTGGTCGGTCGTGGAGGTGGTGATCATCGCGGTGGTGGTGATGATCACCCATTGCTACTACGGGTACACCGCCAGCGGTGGCCCGGTTGGGGTCGGCCAGGCGGTTGGTCGATCGATGCGTTTCTCGCTGGTCTCGGTGGTGGTCGTTGTCCTGCTGGCCGAGTTGGCGCTCTACGGCGTCGACCCGAACTTCAATCTCACGGTGTAGCCGCGGTGCCAACGCTGGTGACGAGGAAGAACCGACGTGCGTGGCTGTATGTGGAGGGTGTTGTCCTGCTGTTGGTGGGCGCGTTGGTGCTCGTATTGGTGTACAAGCAGTTTCGTGGGGAATTCACGCCGAAGACCGAGCTGACTATGGTCGCCTCCCGGGCTGGGCTGGTTATGGAAGCTGGATCCAAAGTCACCTACAACGGGGTGGAGATCGGCCGGGTGGGCAGCATTTCGGAGATTGAGCGTGACGGCCGGCCGGCGGCGAAGCTGGTTTTGGACGTGAATCCTCGCTACATCAGCCTGATTCCGGTCAATGTGGTGGCCGATATCGAGGCGGCCACCCTGTTCGGCAACAAGTATGTTGCGCTGTCCGCGCCGAAAATTCCTCAACAGCAGCGGATTTCCTCACATGACGTGATTGATGTGGGGTCGGTGACCACCGAATTCAACACGTTGTTCGAGACGATCACCTCGATCGCCGAGAAGGTGGATCCGATCGAGCTGAACGCGACGCTGTCCGCGGTAGCACAGGCGCCGGATGGGCTGGGCGGCAAGTTCGGTGAGTCGATCGTTAATGGCAATCAGATTCTGGCGCAATTAAATCCGCGGCTGCCGCAGCTCGGCTATGATGTTCGGCGGTTGGCGGATCTCGGTGAGGTCTATGTCGATGCTTCGCCGGATCTGTGGTCCTTTCTGCAGAACGCACTGACCACTGCGCGCACATTGACCAGCCAACAGCGCGATCTGGATGCCGCGTTGTTGGCGGCTACGGGTGCGGGCAACACCGGTGAAGACGTTTTTGCTCGAGGCGGGCCGTATCTTGCGCGCGCAGCCGCCGATCTGGTGCCCACCGCTACGCTGCTGGACACCTACAGTCCCGAACTGTTCTGCATGATCCGCAACTTTCACGACGCTGCGCCCAAAGTCGCGGACGCGGTGGGCGGCAACGGCTATTCGCTAGCGGCCGCCGGAACGATTTTGGGAGCACCCAATCCCTATGTCTATCCGGACAATCTGCCGCGGGTGAATGCCCACGGTGGACCCGGGGGCCGACCGGGCTGCTGGCAGACGATCACCCGGGAGCTGTGGCCGGCACCCTATCTGGTGATGGACACCGGTGCCAGCCTCGCACCGTACAACCACGTCGAGCTCGGCCAACCGATGTTCACTGAATACGTATGGGGACGCCAATACGGAGAGAACACGATCAACCCATGAAAACCACAGGCACAACTATCAAACTCGGCATCGTCTGGTTGGTGCTGTCGGTGTTCACCGTGATGATCATCGTGGTGTTCGGGCAGGTGCGGTTCCATCACACCACCGGGTACTCCGCGGTGTTCACCCATGTCAGCGGGCTGCGGGCCGGGCAATTTGTCCGCGCTGCGGGCGTAGAGGTCGGCAAGGTCGCCAAGGTAACGCTGATCGACGGGGACAAGCAAGTATTGGTGGACTTCACCGTGGATCGCTCGCTGTCACTGGATCAGGCGACGACCGCCTCGATCCGCTACCTCAACCTGATCGGCGACCGGTACCTTGAGCTCGGCCGCGGTCACAGCGGTCAGCGGCTGGCGCCGGGTGCCACGATCCCGCTCGAGCACACCCATCCGGCCTTGGATCTCGACGCTCTGCTCGGCGGGTTTCGCCCACTCTTCCAAACGTTGGACCCAGACAAGGTCAACAGCATCGCCTCCTCGATCATCACCGTGTTCCAAGGGCAAGGCGCCACCATCAACGACATCCTCGACCAGACCGCCTCGCTGACGGCAACGCTGGCCGACCGGGACCATGCGATAGGTGAGGTCGTCAACAACTTGAACACCGTGCTGGCCACCACCGTCAAGCATCAAACGGAATTCGACCGCACGGTCGACAAGCTAGAGGTGCTGATCACTGGACTGAAGAACAGGGCGGACCCGCTGGCCGCGGCGGCGGCACACATCAGCAGCGCCGCGGGAACCCTAGCCGACCTGCTGGGGGCGGATCGTCCATTGCTGCACAGCAGCTTCGGGCACCTCGAGGGCATCCAGCAGCCGCTCATAGACGAGCTGGCAGAACTCGACCACGTGTTGGGCAAGCTGCCGGACGCCTACCGGATCATCGGCCGCGCCGGCGGCATATACGGTGACTTCTTCAACTTCTATCTGTGTGACATCTCACTGAAAGTCAACGGATTACAGCCTGGAGGTCCGGTACGCACCGTCAAGTTGTTCGGCCAGCCGACCGGCAGGTGCACACCGCAATGAGAACGCTGACCGAGTTCAACCGCGGCCGTGTCGGGATGATGGGTGCGGTGGTCACGGTGCTCGTCGTTGGTGTTGCGCAAAGCTTCACCAGCGTGCCGATGCTGTTCGCCACACCTACCTACTATGCGCAATTCGCCGACATGGGTGGCATCAACACGGGCGATAAGGTGGAAATCGCTGGGGTGAACGTCGGGCTGGTGCGCTCGCTGGCAATCCGCGGCAACCGCGTGTTGATCGGATTCTCGTTGCCCGGCAAGACAATCGGGATGCAAAGCCGGGCAGCAATTCGCACCGACACCATTCTTGGCCGTAAGAACCTGGAGATCGAACCCCGCGGTTCGGAGCCGTTGAAACCCAACGGTTTCCTGCCGTTGGCGCAGACCACTACGCCATACCAAATCTATGACGCGTTCGTCGATGTCACGAAGGCGGCGACGGGCTGGGACATCGATGCCGTCAAACGCTCGCTAAACGTGTTGTCGGAGACATTCGATCAGACCGCCCCGCATCTAAGTGCCGCCCTCGAGGGTGTCAAGGCATTCTCCGACACCGTCGGCCGGCGCGGCGAGCAGATCGAGCAACTGCTGGCGAACGCCAACAGGATCGCGCGCGTGCTCGGCGACCGCAGCGAGCAGGTCAACGGGCTGCTGGTGAATGCCAAGACGCTGCTGGCCGCGTTCAAGCAACGCAGCCAGGCACTGCGCATTCTGCTAACCAACGTGTCGGAGGCATCAGCCCAGGTATCTGGCCTGATCACAGACAACCCCAACCTCAACCATGTGCTGGCCCAGTTGCGCACGGTCAGCGAGGAGCTGGTGAAGCGCAAGAACGAATTGGCCGATGTAGCCGTCTTGCTCGGCAGATACACCGCGGCCCTGACAGAGGCCGTCGGTTCCGGACCGTTCTTCAAGGCGATGGTGGTCAATCTGCTGCCCTACCAGATTCTTCAGCCCTGGGTTGACGCGGCGTTCAAAAAGCGGGGCATCGACCCGGAGAACTTCTGGCGCAGTGCGGGTCTGCCGGAATTCCGCTGGCCCGACCCCAACGGCACCCGGTTCCCCAACGGCGCGCCGCCGGCGGCGCCACCGGTGCGGGAGGGTACACCCAAGCATCCGGGACCGGCCGTCCCGCCGGGAACGCCGTGCTCCTACACACCGGCGGCGGGCGCGTTGCCACGGCCCGACAACCCACTACCCTGCGCGGGCGCCACCGTTGGCCCGTTCGGTGGACCCGACTTCCCGGCACCGCTCGATGTCCAGCCGTCGCCGCCTAATCCCGATGGGCCGCCGCCGACGCCGGGCATCCTAAGTGCTGGGCGGCCGGGCGAGCCGGCTCCGGCTGTTCCGGGCATACCGATGCCGCTGCCGCCGAACGCGCCGCCGGGTGCCCGCACCCAACCGCTGGAGCCGTTTCCTGACGGGACGGGAGGTAGCAACCAATGAGCACCATCTTCGACATCCGCAGCCTGCGACTGCCGAAACTGTCTGCAAAGGTAGTGGTCGTCGGCGGGTTGGTGGTGGTCTTGGCGGTCGTGGCCGCTGCGGCCGGCGCGCGGCTCTACCGGAAACTGACTACCACTACCGTGGTCGCGTATTTCTCTGAGGCGCTCGCGCTGTACCCAGGAGACAAAGTCCAGATCATGGGTGTGCGGGTCGGTTCTATCGACAAGATCGAGCCGGCCGGCGACAAGATGCGAGTCACGTTGCACTACAGCAACAAATACCAGGTGCCGGCCACGGCTACCGCGTCGATCCTCAACCCCAGCCTGGTGGCCTCGCGCACCATCCAGCTGTCACCGCCGTACACCGGCGGCCCGGTCTTGCAAGACGGCGCGGTGATCCCAATCGAGCGCACCCAGGTGCCCGTCGAGTGGGATCAGTTGCGCGATTCCATCAATGGGATCCTCCGCCAGCTCGGCCCGACGGAGCGGCAGCCGAAGGGGCCGTTCGGCGACCTCATCGAATCGGCCGCGGACAACCTGGCCGGCAAGGGCAGGCAGCTCAACGAAACGCTGAACAGTTTGTCGCAGGCGTTGACCGCGCTGAACGAGGGCCGGGGAGACTTCGTTGCGATCACGCGAAGCCTGGCGCTATTTGTCAGCGCGCTCTACCAGAATGATCAACAGTTCGTTGCGCTCAACGAAAACCTTGCCGAGTTCACCGACTGGTTCACCAAATCCGACCATGACTTGGCCGACACGGTGGAACGGATCGACGACGTTCTCGGCACCGTCCGAAAGTTCGTGAGCGACAACAGATCCGTGCTGGCTGCCGATGTCAACAACCTCGCCGACGCGACCACTACACTAGTGCAACCCGAGCCGCGGGACGGTCTGGAAACCGCGTTGCACGTGTTGCCGACCTACGCCAGCAACTTCAACAACCTTTACTATCCACTGCACAGCTCTCTGGTGGGCCAGTTCGTGTTCCCCAACTTCGCGAACCCAATTCAGCTCATTTGCAGCGCTATTCAGGCCGGCAGCCGACTCGGCTATCAGGAATCCGCCGAGCTGTGCGCGCAGTACTTGGCACCGGTTCTGGACGCTCTCAAGTTCAATTACTTGCCGTTCGGCTCAAACCCGTTCAGTTCGGCGGCCACTTTGCCCAAGGAGGTGGCTTACTCCGAGGAGCGGCTCCGCCCGC",
    "contigLength": 4349904,
    "proteins": {
     "Rv0589": [
      "VPTLVTRKNRRAWLYVEGVVLLLVGALVLVLVYKQFRGEFTPKTELTMVAFRAGLVMEAGSKVTYNGVEIGRVGSISEIERDGRPAAKLVLDVNPRYISLIPVNVVADIEAATLFGNKYVALSAPKIPQQQRISSHDVIDVGSVTTEFNTLFETITSIAEKVDPIELNATLSAVAQALDGLGGKFGESIVNGNQILAQLNPRLPQLGYDVRRLADLGEVYVDASPDLWSFLQNALTTARTLTSQQRDLDAALLAATGAGNTGEDVFARGGPYLARAAADLVPTATLLDTYSPELFCMIRNFHDAAPKVADAVGGNGYSLAAAGTILGAPNPYVYPDNLPRVNAHGGPGGRPGCWQTITRELWPAPYLVMDTGASLAPYNHVELGQPMFTEYVWGRQYGENTINP*",
      1215
     ],
     "Rv0590": [
      "MKTTGTTIKLGIVWLVLSVFTVMIIVVFGQVRFHHTTGYSAVFTHVSGLRAGQFVRAAGVEVGKVAKVTLIDGDKQVLVDFTVDRSLSLDQATTASIRYLNLIGDRYLELGRGHSGQRLAPGATIPLEHTHPALDLDALLGGFRPLFQTLDPDKVNSIASSIITVFQGQGATINDILDQTASLTATLADRDHAIGEVVNNLNTVLATTVKHQTEFDRTVDKLEVLITGLKNRADPLAAAAAHISSAAGTLADLLGRIVHCCTAASGTSRASSSRS*",
      828
     ],
     "Rv0590A": [
      "LLHSSFGHLEGIQQPLIDELAELDHVLGKLPDAYRIIGRAGGIYGDFFNFYLCDISLKVNGLQPGGPVRTVKLFGQPTGRCTPQ*",
      255
     ],
     "Rv0591": [
      "MRTLTEFNRGRVGMMGAVVTVLVVGVAQSFTSVPMLFATPTYYAQFADTGGINTGDKVEIAGVNVGLVRSLAIRGNRVLIGFSLPGKTIGMQSRAAIRTDTILGRKNLEIEPRGSEPLKPNGFLPLAQTTTPYQIYDAFVDVTKAATGWDIDAVKRSLNVLSETFDQTAPHLSAALEGVKAFSDTVGRRGEQIEQLLANANRIARVLGDRSEQVNGLLVNAKTLLAAFKQRSQALRILLTNVSEASAQVSGLITDNPNLNHVLAQLRTVSEELVKRKNELADVAVLLGRYTAALTEAVGSGPFFKAMVVNLLPYQILQPWVDAAFKKRGIDPENFWRSAGLPEFRWPDPNGTRFPNGAPPAAPPVREGTPKHPGPAVPPGTPCSYTPAAGALPRPDTPLPCAGATVGPFGGPDFPAPLDVQPSPPNPDGPPPTPGILSAGRPGEPAPAVPGIPMPLPPNAPPGARTQPLEPFPDGTGGSNQ*",
      1446
     ]
    }
   }
  },
  "copies": {
   "genome": "E. coli O157:H7 Sakai",
   "base": "Sakai",
   "contig": "Sakai",
   "kind": "copies",
   "prot": "b0497",
   "name": "rhsD",
   "before": [
    {
     "id": "Sakai_00577",
     "rvs": "b0496",
     "s": 614214,
     "e": 616628,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_b0497",
     "rvs": "b0497",
     "s": 617058,
     "e": 621254,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_02074",
     "rvs": "",
     "s": 2043875,
     "e": 2045983,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_b0497",
     "rvs": "b0497",
     "s": 2046050,
     "e": 2050252,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_02076",
     "rvs": "",
     "s": 2050264,
     "e": 2050449,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_02077",
     "rvs": "",
     "s": 2050449,
     "e": 2050658,
     "st": "+",
     "pseudo": false
    }
   ],
   "refResidues": 1426,
   "loci": [
    [
     616500,
     621800
    ],
    [
     2045500,
     2050800
    ]
   ],
   "genes": [
    {
     "id": "Sakai_00577",
     "rvs": "b0496",
     "s": 614214,
     "e": 616628,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_b0497",
     "rvs": "b0497",
     "s": 617058,
     "e": 621254,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_02074",
     "rvs": "",
     "s": 2043875,
     "e": 2045983,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_b0497_2",
     "rvs": "b0497",
     "s": 2046050,
     "e": 2050252,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_02076",
     "rvs": "",
     "s": 2050264,
     "e": 2050449,
     "st": "+",
     "pseudo": false
    },
    {
     "id": "Sakai_02077",
     "rvs": "",
     "s": 2050449,
     "e": 2050658,
     "st": "+",
     "pseudo": false
    }
   ],
   "bakta": [
    {
     "tag": "Sakai_00577",
     "s": 614214,
     "e": 616628,
     "st": "+"
    },
    {
     "tag": "Sakai_00578",
     "s": 617058,
     "e": 621254,
     "st": "+"
    },
    {
     "tag": "Sakai_00579",
     "s": 621235,
     "e": 621495,
     "st": "+"
    },
    {
     "tag": "Sakai_02074",
     "s": 2043875,
     "e": 2045983,
     "st": "+"
    },
    {
     "tag": "Sakai_02075",
     "s": 2046050,
     "e": 2050252,
     "st": "+"
    },
    {
     "tag": "Sakai_02076",
     "s": 2050264,
     "e": 2050449,
     "st": "+"
    },
    {
     "tag": "Sakai_02077",
     "s": 2050449,
     "e": 2050658,
     "st": "+"
    }
   ]
  },
  "copies-1": {
   "genome": "E. coli O157:H7 Sakai",
   "base": "Sakai",
   "contig": "Sakai",
   "from": 616700,
   "to": 621600,
   "bakta": [
    {
     "tag": "Sakai_00578",
     "s": 617058,
     "e": 621254,
     "st": "+",
     "name": "rhsD",
     "product": "protein RhsD",
     "case": 4
    },
    {
     "tag": "Sakai_00579",
     "s": 621235,
     "e": 621495,
     "st": "+",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    }
   ],
   "ref": [
    {
     "prot": "b0497",
     "name": "rhsD",
     "s": 617058,
     "e": 620825,
     "st": "+",
     "len": 1427,
     "match": 3549,
     "num": 3783,
     "complete": false
    }
   ],
   "bad": [
    {
     "prot": "bakta_ref_00534",
     "s": 620544,
     "e": 621251,
     "st": "+"
    },
    {
     "prot": "bakta_ref_00535",
     "s": 621235,
     "e": 621492,
     "st": "+"
    }
   ],
   "refMatch": [
    {
     "tag": "Sakai_00579",
     "match": 255,
     "num": 258
    }
   ],
   "amap": [
    {
     "id": "Sakai_b0497",
     "rvs": "b0497",
     "s": 617058,
     "e": 621254,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "Sakai_b0497",
     "rvs": "b0497",
     "s": 617058,
     "e": 620825,
     "st": "+"
    }
   ],
   "translate": {
    "offset": 615499,
    "dna": "TCGCCGCATGACGCTGAAATCGCTGCCTCTGCGCCTGGCGGTTAGCCGCCTGTTACGTCAGCCGTGGTCAACGTTAAGCCAGCTTTCGGCATTTTCGCTCTCCTTTATGCTGCTGGCACTGCTACTGGTGCTGCGTGGCGATCTGCTCGACCGCTGGCAGCAGCAGTTACCGCCAGAAAGCCCGAACTACTTTTTAATCAACATCGCCACAGAACAGGTTACGCCGCTAAAAGCGTTCCTCGCGGAACATCAGATAGTCCCGGAATCGTTTTATCCGGTGGTGCGAGCGCGGCTGACGGCGATTAACGGCAAGCCGACAGAAGGCAATGAAGATGAGGCGCTCAACCGCGAACTCAATCTTACCTGGCAAAATACGCGGCCCGATCATAATCCGATTGTCGCCGGTAACTGGCCGCCAAAAGCTGATGAAGTTTCGATGGAAGAGGGGCTGGCGAAACGCTTAAACGTTACCCTCGGCGATACCGTGACTTTTATGGGCGATACCCAGGAGTTCCGCGCTAAAGTGACCAGCCTGCGCAAAGTGGACTGGGAAAGTCTGCGGCCTAATTTCTATTTTATTTTCCCTGAAGGGGCATTAGACGGGCAACCGCAGAGCTGGCTTACCAGTTTCCGCTGGGAGAATGGCAACGGCATGTTGACGCAACTCAACCGCCAGTTCCCTACCATTAGCCTGTTAGACATTGGCGCGATTTTAAAACAGGTCGGTCAGGTGCTGGAGCAGGTAAGTCGGGCGCTGGAGGTGATGGTGGTGCTGGTCACCGCCTGCGGTATGTTGCTGTTGCTGGCACAGGTGCAGGTGGGAATGCGTCAGCGTCATCAGGAGCTGGTGGTGTGGCGCACACTCGGTGCGGGGAAAAAACTGCTGCGTACCACGTTGTGGTGTGAGTTCGCCATGCTTGGGTTTGTTTCCGGCCTGGTGGCCGCAATTGGTGCGGAAACGGCACTGGCAGTGTTGCAGGCGAAAGTGTTTGATTTTCCGTGGGAGCCAGACTGGCGATTGTGGATTGTTCTGCCGTGCAGCGGTGCGCTGCTGCTGTCGCTTTTCGGCGGCTGGCTGGGCGCGAGACTGGTTAAGGGTAAGGCGCTGTTCAGGCAGTTTGCGGGGTGATGAGAGTAAATAAGGGGCGTTCTGGTTGAATGGAACGCCTTGTTTAAATTAATTCTACAATATGTTTCATTGATGATGTTATTGAATTGGTGATTTCCTATCTTTCTATTGCTGATATTAATCTGAATCGGTGGGGTTTAATAATTTTAAATATAGCCCTACAAAAATAACTCAGACTAATAAAATCATAAATCATATGCGTTGAATGGATATTATCCATATAGTGAATTTGTTGATGATGAATTCATCTGTGCTAAAAATTTTAGTTTAATAAAATATTGAAAGTGACCTGTAATAACAGTTGTTGTTGATTGAGAACAAATAAGTTTATGTGAAAAATATATAAATACATTAGCTGGTCTTGTGTGTCATTTTATTTTTTTTGTTGCTAACACAGGGATATGAACAATAACTAAAAGGGCACTTTATATGAGCGGAAAACCAGCGGCGCGTCAGGGAGATATGACTCAGTATGGCGGTCCCATTGTCCAGGGTTCGGCAGGTGTAAGAATTGGCGCGCCAACTGGCGTGGCGTGCTCGGTGTGCCCGGGCGGGATGACTTCGGGCAACCCGGTAAATCCGCTGCTGGGGGCGAAGGTGCTGCCCGGCGAGACGGACCTTGCGCTGCCCGGCCCGCTGCCGTTCATTCTCTCCCGCACCTACAGCAGCTACCGGACCCGGACGCCTGCGCCGGTGGGGATTTTTGGCCCCGGCTGGAAAGCGCCTTCTGATATCCGCTTACAGCTACGCGATGATGCACTGGTACTCAATGACAACGGCGGGCGGAGCATTCACTTTGAGCCGCTGCTGCCGGGGGAGGCGGTGTACAGCCGCAGCGAGTCAATGTGGCTGGTGCGCGGTGGTAAGGCAGCGCAGCCGGACGGCCACACGCTGGCGCGGCTGTGGGGGGCGCTGCCGCCGGATATCCGGTTAAGCCCGCATCTTTACCTGGCGACCAACAGCGCACAGGGGCCGTGGTGGATACTGGGGTGGTCAGAGCGGGTGCCGGGTGCTGAGGACGTACTGCCAGCGCCGCTGCCGCCGTACCGGGTGCTTACCGGGATGGCGGACCGCTTCGGGCGGACGCTGACGTACAGGCGTGAGGCCGCCGGTGACCTGGCCGGGGAAATCACCGGCGTGACGGACGGTGCCGGGCGGGAGTTCCGTCTGGTGCTGACCACGCAGGCGCAGCGGGCGGAAGAGGCCCGTAAACAGCACACCGCTTCTTTATCTTCCCCTGACACCCCCCGCCCTCTTTCAGACTCAGCGTTCCCCGACACACTGCCCGGTACCGAATACGGTCCCGACAGAGGTATCCGCCTTTCGGCGGTGTGGCTGACGCACGACCCGGCATACCCGGAGAGCCTGCCCGCTGCGCCACTGGTGCGGTACACGTATACGGAAGCCGGTGAACTGCTGGCGGTATATGACCGCAGCAATACGCAGGTGCGCGCTTTCACGTATGACGCGCAGCACCCGGGCCGGATGGTGGGGCACCGTTATGCGGGAAGGCCGGAGATGCGCTACCGCTACGACGATGCGGGGCGGGTGGTGGAGCAACTGAACCCGGCAGGCCTGAGTTACCACTACCAGTATGAGCAGGACCGCATCACCGTCACGGACAGCCTGAACCGGCGTGAGGTGCTGCATACAGAAGGCGGGGCCGGGCTGAAGCGGGTGGTGAAAAAAGAACTGGCGGACGGCAGCGTCACGCACAGCGGGTATGACGCGGCAGGAAGGCTAACGGCGCAGACGGACGCGGCGGGACGGCGGACAGAGTATGGTCTGAATGTGGTATCCGGCGATATCACGGACATCACCACACCGGACGGGCGGGAGACGAAATTTTACTATAACGACGGGAACCAGCTGACGGCGGTGGTGTCCCCGGACGGGCTGGAGAGCCGCCGGGCATATGATGAACCGGGCAGGCTGGTATCGGAGACATCGCGCAGCGGGGAGACAGTACGCTACCGCTACGATGATGCGTACAGTGAGTTACCGGCGACGACAACAGATGCGACGGGCAGCACCCGGCAGATGACCTGGAGCCGCTACGGTCAGTTGCTGGCGTTCACCGACTGCTCGGGCTACCAGACCCGCTATGAATACGACCGCTTCGGCCAGATGACGGCGGTCCACCGCGAGGAAGGCATCAGCCTTTACCGCCGCTATGACAACCGTGGCCGGTTAATCTCGGTGAAAGACGCACAGGGCCATGAAACGCGGTATGAGTACAACGCCGCAGGCGACCTGACTGCCGTTATCACCCCGGACGGCAACCGGAGCGAGACACAGTACGATGCGTGGGGAAAGGCGGTCAGCACCACGCAGGGCGGGCTGACGCGCAGTATGGAATACGACCTTGCCGGACGCATCACCACGCTGACCAACGAGAACGGCAGCCGGAGTGAGTTTACCTACGATGCGCTTGACCGGCTGGTACAGCAGCGCGGCTTTGACGGGCGGACGCAACGTTACCACTATGACCTGACCGGAAAACTCACGCAGAGTGAAGATGAGGGGCTTGTCACCCTCTGGCACTACGACGAATCGGACCGCCTCACTCACCGCACGGTGAACGGCGAACCGGCAGAGCAGTGGCAGTACGACGAGCACGGCTGGCTGACAGAAATCAGCCACCTGAGCGAAGGCCATCAGGTGGCGGTGCATTACGGTTATGATGATAAGGGCCGCCTGGCCGGGGAGCGCCAGACGGTGCATAACCCGGAGACGGGGGAACTGCTGTGGCAGCATGAGACAGAGCACGCATACAACGAACAGGGTCTGGCAAACCGCGTCACGCCGGACAGCCTGCCGCGGGTGGAGTGGCTGACCTACGGCAGCGGTTATCTTGCGGGGATGAAGCTGGGCGGGACGCCGCTGGTGGAGTTCACGCGCGACAGGCTGCACCGCGAGACGGTGCGCAGCTTCGGCAATAACGCATACGAACTGACCAGCACATACACTCCCGCAGGCCATTTACAGAGCCAGCGCCTGAACAGCCAGGTGTATGACCGTGACTACGACTGGAATGACAATGGCGACCTGGTGCGCATCAGCGGCCCGCGACAGACGCGGGAATATGGCTACAGCGCCACGGGCAGGCTGGAGAGCGTGCGCACCCTTGCATCAGACCTGGATATCCGCATCCCGTATGCGACCGACCCGGCGGGAAACCGGCTGCCGGACCCGGAGCTACACCCGGACAGCACGCTCACGGCGTGGCCGGATAACCGCATCGCGGAGGATGCGCACTATGTCTACCGACACGATGAATACGGCAGGCTGACGGAGAAGACGGACCGCATCCCGGCGGGTGTGATACGGACGGACGACGAGCGGACCCACCACTACCACTACGACAGCCAGCACCGCCTGGTGTTCTACACGCGGATACAGCATGGCGAGCCACTGGTCGAGAGCCGCTACCTCTACGACCCGCTGGGCCGCCGGACGGGGAAACGGGTGTGGCGTCGCGGGCGTGACCTGACGGGATGGATGTCGCTGTCGCGGAAACCGGAGGTGACGTGGTACGGGTGGGACGGCGACCGGCTGACGACGGTACAGACCGACACCACGCGTATCCAGACGGTATACGAGCCGGGAAGCTTCGCGCCGCTCATCCGCATCGAAACAGACAACGGTGAGCGGGAGAAAGCGCAGCGCCGCAGCCTGGCGGAGAAGCTGCAGCAGGAAGGGAGCGAGGACGGGCACGGTGTGGTATTTCCGGCTGAACTGGTGCGGCTGCTGGACAGACTGGAGGAAGAAATCCGGGCAGACCGCGTGAGCAGTGAAAGCCGGGCGTGGCTTGCGCAGTGCGGGCTGACGGTGGAGCAACTGGCCAGACAGGTGGAGCCGGAATACACACCGGCGAGAAAAGTTCATCTTTACCACTGCGACCACCGGGGCCTGCCGCTGGCGCTCATCAGCGAAGACGGCAATACGGCGTGGAGCGGGGAGTATGATGAATGGGGCAACCAGCTGAATGAGGAGAACCCGCATCACCTGCACCAGCCGTACCGGCTGCCGGGGCAGCAGTATGATAAGGAGTCGGGGCTGTACTACAACCGGAACCGGTACTACGATCCGTTGCAGGGGCGGTATATCACTCAGGACCCGATAGGGCTGGAGGGGGGATGGAGTCTGTATGCGTATCCGCTGAATCCGGTGAATGGTATTGATCCATTAGGGTTAAGTCCCGCAGATGTAGCGCTAATAAGAAGAAAAGATCAACTAAACCATCAAAGAGCATGGGATATATTATCTGATACTTATGAAGATATGAAGAGATTAAATTTAGGTGGGACTGATCAATTTTTCCATTGTATGGCATTTTGTCGAGTGTCTAAATTAAATGACGCTGGTGTTAGCCGATCGGCGAAAGGGCTGGGTTATGAAAAAGAGATTAGAGATTACGGGTTAAATCTGTTCGGTATGTACGGCAGAAAAGTAAAGCTATCCCATTCTGAAATGATTGAAGATAATAAAAAAGACTTGGCTGTAAATGACCATGGGTTGACATGTCCATCAACAACAGATTGCTCAGATAGATGTAGTGATTATATTAATCCAGAGCATAAAGAAACGATAAAGGCTTTACAAGATGCTGGCTATCTCAAGTAATCTATCAAAGATGATAATATTTATTATTGCTATTATAATCATTGTTGTTTTATGCGTAATTACTTATCTTTATTTATACAAAGATGAATCTCTTGTAAGTAAACATTACATAAACTATATGGCAATACCAGAAAATGATGGAGTTTTTACATGGCTCCCAGATTTTTTTCCGCACGTAGCGGTGGATATATCAATATACACAAATGTAGAAGATGATTATTTTTTTCTTATTTTTCCCTAACAAATGATGATGGGGGGAGGTTTAAGAAAACATTGACAGTGAGGGCCAGGGAACAAGTGGCGAAAATCGTATCAAAGAATGATCCAGATACAAAAAAAGTGTGGTGTAAATATGGTAAGATACCAGGGCAAGGGGATGGTGTAAACCTTTTTTTGTTGGTGAAATTAATGTTACGCATTATTTTATAACAAATATTGGAGCTGGATTGCCTGATGCTTGTGCAGAGTAATTGCTTGAATTAAGAGTCTATCCCATATCGAAGTCGTCAACTTCGTAGTGAGGAAAAGTAAAATTCCTGACTGAGAAAAGACATGTCGGCTATTGTGTAAAGCCATATAGCTCAGACGATGAATATCTACTCGTATTCAGTTGTTTATTGAGGGTGAGTTCCGACCCTGAAACAACAAATAAAATGAACAGTCAGAGAGTTTACATAGAATTGCACTGGTCTTTTACGATATCTGACATTGTGTAATACATATTCAGCCATGCATTAATTAAAATGTTACGTGTTTAATGTGAGTCCCTATCTGAAAATAAATAATCCTTCCGGATTAAAATAAATTCTTGCCGGGAAAGAAAGAGGAAATAAACCATTAGCGGAAAACCAGCGGCACGCCAGGGTGACATGCCCCAGTACGGCGGCCCGATAGTCCAGGGTTCAGCGGGGGTACTGATAGTCGCGCCGACTGCATCGAGAACTCCATACAACTTTTTACGGAACGTCACCTGATAAAATTTGTTCAGTATCGTCTTTTGAGATTCTCAAGGTTATGGAGTAACCAGACAGAACGGACACCACTGTCGGAGATAAAAACCCCCTGTTTACGCAGTTCGTTGCTGGCCCGGTGCTGACCATGTGTCGGGAAAGCGACGGCGTAATCAACAACAGCCTGTTCAGTTGCCTCATCGGTACGGTTCTTAAGGTTAGGAGCGCGGCGACTACGATTAATCTGCGCATCCACACCGCCTTCAGCGACCAGTTCGCGGTAACGATAAAACGTATCACGCGAAACGCCCATGATTTTACAGGCTTTTGATACGTTGCTGAGTTCTTCAGCCAGATTGAGCAAACCGGCTTTGTGTTTGATGACGGGATTGGCAGTATGAAGCATGAGAGTTACCTCTTGTTTTGGATAAGGATTCGACACTCATATCAAAACCGGTAACTCTCAACCTTTCAAGGCCATGTGTCAGATCAAGTCGCGACTAATACAAATACGTCCCTCATTACTGCGCCTTAACCCATTCCGCCACTTCCGCCCACTCACCGCGAAAGACAACTTTTTCCGCTTTTTTCTCAAG",
    "contigLength": 5498578,
    "proteins": {
     "b0497": [
      "MSGKPAARQGDMTQYGGPIVQGSAGVRIGAPTGVACSVCPGGMTSGNPVNPLLGAKVLPGETDLALPGPLPFILSRTYSSYRTKTPAPVGVFGPGWKAPSDIRLQLRDDGLILNDNGGRSIHFEPLLPGEAVYSRSESMWLVRGGKAAQPDGHTLARLWGALPPDIRLSPHLYLATNSAQGPWWILGWSERVPGAEDVLPAPLPPYRVLTGMADRFGRTLTYRREAAGDLAGEITGVTDGAGREFRLVLTTQAQRAEEARTSSLSSSDSSRPLSASAFPDTLPGTEYGPDRGIRLSAVWLMHDPAYPESLPAAPLVRYTYTEAGELLAVYDRSNTQVRAFTYDAQHPGRMVAHRYAGRPEMRYRYDDTGRVVEQLNPAGLSYRYLYEQDRITVTDSLNRREVLHTEGGAGLKRVVKKELADGSVTRSGYDAAGRLTAQTDAAGRRTEYGLNVVSGDITDITTPDGRETKFYYNDGNQLTAVVSPDGLESRREYDEPGRLVSETSRSGETVRYRYDDAHSELPATTTDATGSTRQMTWSRYGQLLAFTDCSGYQTRYEYDRFGQMTAVHREEGISLYRRYDNRGRLTSVKDAQGRETRYEYNAAGDLTAVITPDGNRSETQYDAWGKAVSTTQGGLTRSMEYDAAGRVISLTNENGSHSVFSYDALDRLVQQGGFDGRTQRYHYDLTGKLTQSEDEGLVILWYYDESDRITHRTVNGEPAEQWQYDGHGWLTDISHLSEGHRVAVHYGYDDKGRLTGECQTVENPETGELLWQHETKHAYNEQGLANRVTPDSLPPVEWLTYGSGYLAGMKLGGTPLVEYTRDRLHRETVRSFGSMAGSNAAYELTSTYTPAGQLQSQHLNSLVYDRDYGWSDNGDLVRISGPRQTREYGYSATGRLESVRTLAPDLDIRIPYATDPAGNRLPDPELHPDSTLTVWPDNRIAEDAHYVYRHDEYGRLTEKTDRIPAGVIRTDDERTHHYHYDSQHRLVFYTRIQHGEPLVESRYLYDPLGRRMAKRVWRRERDLTGWMSLSRKPEVTWYGWDGDRLTTVQTDTTRIQTVYEPGSFTPLIRVETENGEREKAQRRSLAETLQQEGSENGHGVVFPAELVRLLDRLEEEIRADRVSSESRAWLAQCGLTVEQLARQVEPEYTPARKAHLYHCDHRGLPLALISEDGNTAWSAEYDEWGNQLNEENPHHVYQPYRLPGQQHDEESGLYYNRHRYYDPLQGRYITQDPMGLKGGWNLYQYPLNPLQQIDPMGLLQTWDDARSGACTGGVCGVLSRIIGPSKFDSTADAALDALKETQNRSLCNDMEYSGIVCKDTNGKYFASKAETDNLRKESYPLKRKCPTGTDRVAAYHTHGADSHGDYVDEFFSSSDKNLVRSKDNNLEAFYLATPDGRFEALNNKGEYIFIRNSVPGLSSVCIPYHD*",
      4281
     ]
    }
   }
  },
  "copies-2": {
   "genome": "E. coli O157:H7 Sakai",
   "base": "Sakai",
   "contig": "Sakai",
   "from": 2045700,
   "to": 2050600,
   "bakta": [
    {
     "tag": "Sakai_02074",
     "s": 2043875,
     "e": 2045983,
     "st": "+",
     "name": "vgrG",
     "product": "type VI secretion system effector VgrG",
     "case": 1
    },
    {
     "tag": "Sakai_02075",
     "s": 2046050,
     "e": 2050252,
     "st": "+",
     "name": "rhsD",
     "product": "protein RhsD",
     "case": 4
    },
    {
     "tag": "Sakai_02076",
     "s": 2050264,
     "e": 2050449,
     "st": "+",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    },
    {
     "tag": "Sakai_02077",
     "s": 2050449,
     "e": 2050658,
     "st": "+",
     "name": "",
     "product": "hypothetical protein",
     "case": 1
    }
   ],
   "ref": [
    {
     "prot": "b0497",
     "name": "rhsD",
     "s": 2046050,
     "e": 2049892,
     "st": "+",
     "len": 1427,
     "match": 3561,
     "num": 3858,
     "complete": false
    }
   ],
   "bad": [],
   "refMatch": [
    {
     "tag": "Sakai_02074",
     "match": 0,
     "num": 0
    },
    {
     "tag": "Sakai_02076",
     "match": 0,
     "num": 0
    },
    {
     "tag": "Sakai_02077",
     "match": 0,
     "num": 0
    }
   ],
   "amap": [
    {
     "id": "Sakai_02074",
     "rvs": "",
     "s": 2043875,
     "e": 2045983,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "Sakai_b0497_2",
     "rvs": "b0497",
     "s": 2046050,
     "e": 2050252,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "Sakai_02076",
     "rvs": "",
     "s": 2050264,
     "e": 2050449,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    },
    {
     "id": "Sakai_02077",
     "rvs": "",
     "s": 2050449,
     "e": 2050658,
     "st": "+",
     "pseudo": false,
     "partial": "",
     "note": ""
    }
   ],
   "corrected": [
    {
     "id": "Sakai_02074",
     "rvs": "",
     "s": 2043875,
     "e": 2045983,
     "st": "+"
    },
    {
     "id": "Sakai_b0497",
     "rvs": "b0497",
     "s": 2046050,
     "e": 2049892,
     "st": "+"
    },
    {
     "id": "Sakai_02076",
     "rvs": "",
     "s": 2050264,
     "e": 2050449,
     "st": "+"
    },
    {
     "id": "Sakai_02077",
     "rvs": "",
     "s": 2050449,
     "e": 2050658,
     "st": "+"
    }
   ],
   "translate": {
    "offset": 2042674,
    "dna": "AGGATAATGCCACGAATATCAGTGCGATTGCTGATGCAGTATGTCAGTTACCAACGTTACAAGAAATACTAAAAAGGAATGAAATTATTTGAGTGGGATACCAGAAGATTATGTTATTTATTATTCAGTAACTTTATCCACAATGATGTGCGTAATTAATTAGATTCATGTTATGTTTTTTAATTTATATCACACTCCCTTTCATTCAGCTTGTCTATTTCATTTCCTCTGAAAAGAAACTAACTATTCACATTGTGCAGATGTTTCACTTGTTATCACAGGTATTCTATAATCTCAAAATGTTTTTAATGATGGATATGCTCGGAGTTGGAGATGCAATTAATATTAATACAAATAAAAATATCCGGCAGGTATGCTAAAAGGATAAATGGCGTTCTCTGTATTTATTCGCAATTATCTGTAGGACAAATATGGAGAGCTTGATGAGTCAGAGCTAAAACAGCAATTGATTAATAAACATGTTGATTTATTCAAGAAAATGGATGTATTTTAAATATTGATTCTTAATATAACCTAATCAATAGATTATCACGTGAGAGCTGCACATTATTGGAGGGGCAGGGGAAAACACCATTTTGGGAACAATGAATCAGGTACTCTTACCACTGAAGCAAGATCATCCGGGCACCAGGTGGTGGTGCAGTGCTTTGGAAATGGTCCCCTCGGGGGAGTTGTCATTATTTTCTATGCTATTTGCTTATGATATTCATCAGGTGAATTTAATAAATTCCTAATTTGAAATATCCACTATTAAGCCAGTGTTTAACGACGTTAAGATGGCAATGTGCAGATCATGTTATCAGGGTGGGGCGGTAGCGTATATTCCTTATACAAATCATCATTCATTTCAATGATTGTTATTAAATTATTTATTGTTGTGAATAAAAAGACATTATGAATATTGGTTTGACTTTGCATATCAGGTTTTTTATTACTCGCTGTAATGTACGAGCCTGTCGTGAATAACGTCGAGAACCCATCCGTTGCGTTTGTGATTATTTTGTTGACTAAACAGACACCCGTTTCTCTGAAGTAAAATCCCAGACTAAATCATCACATAACCATGACATTTTTCTGATATTCCCCGGTAACGCCAGATGTCGACTCACTTATTTCATCCGCCCGACATGTTGTACATGCAATATGTGCGGATTATTATTCATTTTCACGGAGGTTGCTATGTCAACCGGATTACGTTTCACGCTGGAAGTGGACGGCCTGCCACCGGACGCGTTTGCGGTGGTCTCCTTTCATCTGAACCAGTCACTCTCTTCACTTTTTTCCCTCGACCTTTCTCTGGTCAGTCAGCAGTTTCTTTCCCTTGAATTTGCGCAGGTGCTGGACAAAATGGCCTACCTGACGGTATGGCAGGGCGATGACGTACAGCGCCGGGTGAAAGGTGTGGTGACCTGGTTTGAACTGGGGGAGAACGACAAAAACCAGATGCTGTACAGCATGAAGGTGTGCCCACCGCTGTGGCGCACAGGGCTGCGCCAGAACTTCCGTATCTTCCAGAATGAGGACATCGAAAGCATACTCGCTACGATCCTGAAAGAAAACGGTGTGACCGAGTGGAGTCCGCTGTTCAGCGAGCCGCATCCTTCCCGTGAGTTTTGTGTCCAGTACGGCGAAACTGATTACGATTTCCTGTGCCGGATGGCGGCGGAGGAAGGCATCTTCTTTTATGAGGAGCACGCGCAAAAAAGTATCGACCAGAGCCTGGTCCTGTGCGACACCGTGCGTTATCTGCCGGAGTCCTTTGAGATCCCCTGGAACCCGAACACCCGTACCGAGGTGAGCACCCTATGCATCAGCCAGTTCCGCTACAGCGCACAAATCCGCCCTTCTTCCGTGGTGACCAAGGACTACACCTTTAAACGACCCGGCTGGGCAGGGCGTTTTGATCAGGAAGGCCAGTACCAGGATTACCAGCGCACACAGTATGAAGTGTATGACTACCCCGGACGTTTCAAGGGTGCCCACGGGCAGAACTTTGCCCGCTGGCAGATGGATGGCTGGCGCAACAACGCAGAAGTGGCGCGCGGAACAAGCCGTTCTCCGGAAATATGGCCGGGACGGCGAATTGTGCTGACGGGGCATCCGCAGGCGAACCTGAACCGGGAATGGCAGGTGGTGGCAAGTGAACTGCACGGCGAGCAGCCACAGGCGGTGCCGGGACGCAGGGGTTCAGGTACCACGCTGAATAACCACTTTGCGGTAATACCGGCAGACCGGACATGGCGACCACAGCCGTTGCTGAAACCGCTGGTGGATGGCCCGCAGAGCGCCGTCGTGACGGGACCGGCAGGCGAGGAAATCTTCTGCGACGAACATGGTCGCGTGCGGGTGAAATTTAACTGGGACCGCTATAACCCGTCAAACCAGGACAGTTCATGCTGGATCCGTGTGGCACAGGCGTGGGCAGGCACCGGCTTTGGCAACCTGGCGATACCGCGTGTGGGTCAGGAGGTGATTGTGGACTTCCTCAACGGCGATCCGGACCAGCCGATTATTATGGGGCGTACCTACCACCAGGAAAACCGCACACCCGGCAGCCTGCCGGGGACGAAGACGCAGATGACCATTCGTTCGAAAAACTATAAGGGCAGCGGGTTTAATGAACTGAAGTTTGACGATGCGACAGGGAAAGAACAGGTCTACATCCACGCGCAGAAGAACATGAACACCGAGGTGCTGAATAACCGCACCACGGATGTGATAAACAACCATGCTGAGACCATTGGCAACAATCAGATGATTGCGGTTACCAACAATCAGATACAGACGGTGGGCGTTAACCAGATAGAGACGGTGGGCAGTAACCAGATCATTAACGTGGGTTCTGTTCAGGTGGAAACGATTGGACTTGTTCGTGCGCTGACCGTGGGCGTGGCGTACCAGACGACGGTAGGTGGCATTATGAACACTTCGGTGGCTCTGATGCAGTCCTCGCAGATCGGTTTGCATAAATCGTTGAGGGTCGGGCTGGGTTATGACGTCAAAGTCGGAAATAACGTTACCTTCACCGTTGGTAAAACGAAAAAGGATGATACCGGGCAGACCGCGATTTACTCCGCCGGTGAGCATCTGGAGCTCTGCTGTGGTAAGGCAAGGCTGGTGCTGACGAAGGACGGACAAATTTTTCTCAACGGCACAAAAATTCATTTGCAGGGTAAGGAGCAGGTTAATGGTGACTCACTGTTGATTAACTGGAACTGTGCAGCCTCGAAATCCCCACCGAAGACCCCTGATGAAAAGCAGGATACGCCGGATATGAGAGAGTACTGATCCCGTACCTGTTCGGTATTCAGTAACAAACAAGCATTAAAAATAAACAGAAGAAGGAATTGAATCATGAGCGGAAAACCAGCGGCGCGTCAGGGAGATATGACTCAGTATGGCGGTCCCATTGTCCAGGGTTCGGCAGGTGTAAGAATTGGCGCGCCAACTGGCGTGGCGTGCTCGGTGTGCCCGGGCGGGATGACTTCGGGCAACCCGGTAAATCCGCTGCTGGGGGCGAAGGTGCTGCCCGGCGAGACGGACCTTGCGCTGCCCGGCCCGCTGCCGTTCATTCTCTCCCGCACCTACAGCAGCTACCGGACCCGGACGCCTGCGCCGGTGGGGATTTTTGGCCCCGGCTGGAAAGCGCCTTCTGATATCCGCTTACAGCTACGCGATGATGCACTGGTACTCAATGACAACGGCGGGCGGAGCATTCACTTTGAGCCGCTGCTGCCGGGGGAGGCGGTGTACAGCCGCAGCGAGTCAATGTGGCTGGTGCGCGGTGGTAAGGCAGCGCAGCCGGACGGCCACACGCTGGCGCGGCTGTGGGGGGCGCTGCCGCCGGATATCCGGTTAAGCCCGCATCTTTACCTGGCGACCAACAGCGCACAGGGGCCGTGGTGGATACTGGGGTGGTCAGAGCGGGTGCCGGGTGCTGAGGACGTACTGCCAGCGCCGCTGCCGCCGTACCGGGTGCTTACCGGGATGGCGGACCGCTTCGGGCGGACGCTGACGTACAGGCGTGAGGCCGCCGGTGACCTGGCCGGGGAAATCACCGGCGTGACGGACGGTGCCGGGCGGGAGTTCCGTCTGGTGCTGACCACGCAGGCGCAGCGGGCGGAAGAGGCCCGTAAACAGCACACCGCTTCTTTATCTTCCCCTGACACCCCCCGCCCTCTTTCAGACTCAGCGTTCCCCGACACACTGCCCGGTACCGAATACGGTCCCGACAGAGGTATCCGCCTTTCGGCGGTGTGGCTGACGCACGACCCGGCATACCCGGAAAGCCTGCCCGGTGCGCCACTGGCGCGGTACACGTATACGGAAGCCGGTGAACTGCTGGCGGTATATGACCGCAGCAATACGCAGGTGCGCGCTTTCACGTATGACGCGCAGCATCCGGGCCGGATGGTGGCGCACCGTTACGCGGGAAGGCCGGAGATGCGCTACCGCTACGACGATACCGGGCGGGTGGTGGAGCAGCTGAACCCGGCAGGCCTGAGTTACCGCTACCAGTATGAGCAGGACCGCATCACCGTCACGGACAGCCTGAACCGGCGTGAGGTGCTGCATACAGAAGGCGGGGCCGGGCTGAAGCGGGTGGTGAAAAAAGAACTGGCGGACGGCAGCGTCACGCACAGCGGGTATGACGCGGCAGGAAGGCTCACGGCGCAGACGGACGCGGCGGGACGGCGGACAGAGTACGGTCTGAATGTGGTGTCCGGCGATATCACGGACATCACCACACCGGACGGGCGGGAGACGAAATTTTACTATAACGACGGGAACCAGCTGACGGCGGTGGTGTCCCCGGACGGGCTGGAGAGCCGCCGGGCATATGATGAACCGGGCAGGCTGGTATCGGAGACATCGCGCTGTGGGGACGTCATCCGGTATGCTTATGATAATCCGCACAGTGAATTACCGGCCACGACAACAGATGCGACGGGCAGCACCCGGCAGATGACCTGGAGCCGCTACGGGCAGTTGCTGGCGTTCACCGACTGCTCGGGCTACCAGACCCGTTATGAATACGACCGCTTTGGTCAGATGACGGCGGTCCACCGTGAGGAAGGTATCAGCCGTTACCGCCGCTATGACAACCGTGGCCGGTTAACCTCGGTGAAAGACGCACAGGGCCATGAAACGCGGTATGAGTACAACGCCGCAGGCGACCTGACTGCCGTTATCACTCCGGACGGCAACCGGAGCGAGACACAGTACGATGCGTGGGGAAAAGCGGTCAGCACCACGCAGGGCGGGCTGACGCGCAGTATGGAGTATGACCTCGCCGGACGCATCACCACGCTGACCAACGAGAACGGCAGCCGGAGTGAGTTTACCTACGATGCGCTTGACCGGCTGGTACAGCAGCGCGGCTTTGACGGGCGGACGCAACGTTACCACTATGACCTGACCGGAAAACTCACGCAGAGTGAAGATGAGGGGCTTGTCACCCTCTGGCACTACGACGAATCGGACCGCCTCACTCACCGCACGGTGAACGGCGAACCGGCAGAGCAGTGGCAGTACGACGAGCACGGCTGGCTGACAGAAATCAGCCACCTGAGCGAAGGCCATCAGGTGGCGGTGCATTACGGTTATGATGATAAGGGCCGCCTGGCCGGGGAGCGCCAGACGGTGCATAACCCGGAGACGGGGGAACTGCTGTGGCAGCATGAGACAGAGCACGCATACAACGAACAGGGTCTGGCAAACCGCGTCACGCCGGACAGCCTGCCGCGGGTGGAGTGGCTGACCTACGGCAGCGGTTATCTTGCGGGGATGAAGCTGGGCGGGACGCCGCTGGTGGAGTTCACGCGCGACAGGCTGCACCGCGAGACGGTGCGCAGCTTCGGCAATAACGCATACGAACTGACCAGCACATACACTCCCGCAGGCCATTTACAGAGCCAGCGCCTGAACAGCCAGGTGTATGACCGTGACTACGACTGGAATGACAATGGCGACCTGGTGCGCATCAGCGGCCCGCGACAGACGTGGGAATATGGCTACAGTGCCACGGGCAGGCTGGAGAGCGTGCGCACCCTTGCATCAGACCTGGATATCCGCATCCCGTATGCGACCGACCCGGCGGGAAACCGGCTGCCGGACCCGGAGCTACACCCGGACAGCACGCTCACGGCGTGGCCGGATAACCGCATCGCGGAGGATGCGCACTATGTCTACCGACACGATGAATACGGCAGGCTGACGGAGAAGACGGACCGCATCCCGGCGGGTGTGATACGGACGGACGACGAGCGGACCCACCACTACCACTACGACAGCCAGCACCGCCTGGTGTTCTACACGCGGATACAGCATGGCGAGCCACTGGTCGAGAGCCGCTACCTCTACGACCCGCTGGGACGGCGAATGGCAAAACGGGTCTGGCGGCGGGAGCGTGACCTGACGGGGTGGATGTCGCTGTCGCGTAAACCGGAGGTGACGTGGTATGGCTGGGACGGAGACAGGCTGACGACGGTGCAGACTGACACCACACGTATCCAGACGGTATACGAGCCGGGAAGCTTCACGCCGCTCATCCGGGTCGAGACAGAGAACGGCGAGCGGGAAAAAGCGCAGCGGCGCAGCCTGGCAGAGACGCTCCAGCAGGAAGGGAGTGAGAACGGCCACGGCGTGGTGTTCCCGGCTGAACTGGTGCGGCTGCTGGACAGGCTGGAGGAAGAAATCCGGGCAGACCGCGTGAGCAGTGAAAGCCGGGCGTGGCTTGCGCAGTGCGGGCTGACGGTGGAGCAACTGGCCAGACAGGTGGAGCCGGAATACACACCGGCGCGAAAAGTTCATTTTTACCACTGCGACCACCGGGGCCTGCCGCTGGCGCTCATCAGCGAAGACGGCAATACGGCGTGGCGCGGGGAGTATGATGAATGGGGCAACCAGCTTAATGAAGAGAACCCGTATTACCTGCACCAGCCATACCGTCTGCCGGGGCAGCAGCATGATGAGGAATCAGGGCTGTACTATAACCGGAACCGGTACTATGACCCGCTACAGGGGAGGTATATTACACAAGACCCCATTGGGCTGGCGGGGGGATGGAATCTGTATAATTACCCACTGAATCCGATAATAAGGATGGATCCTTTGGGTTTGTATAATTTATATCAATTATTATATGATGTTTGGCATGATGATTCATATGGAACATCATCAATTGATATTACTGGCAGTGGAGATCTAATATCATTAGGTGGTCATGCAGGACTTGGCGTTGCGTTTGCTAAAAAGAAAGGTGAAATGTTATCTGATATTTGTATTTATGCTACAGCATGCGGACATGCAGGAATTGGTGGTGGGATAAATGCGGCTATCACATATTCTGAGACCAAGTCTTTACCTACATCGGGAGTCAGCAATTCAGTAGGTGTAACGGTTGGCGGCGGAGTTGGGGGGCATTTTGCGTATACTTATGTAGTGGATGTTGATAATCCAGAATCATCGACAGAATCTGTTGGTATCGGTGCAGGTGTTGACGCTTCAGTTATGACTCTGGCTTGTAGAACGTGGCAAGAATGCTGGGTCAATTAACAGGAGGAAATATGAAAAAATATATATTTTTAGTTTTATTTTTATTGTATGGTTGTTCTCCAGTGGACATAGAGCCAGCTAGACAAGAAGTTAATATCTTTCGTGGTTTATATCAGACAGGGAAATACAGTGATATATATAAAATCACTTCAAGTGATTTTCAGACTGCCACATCTGAAAAAAATTTATCGACATAATGAACGAAGCAAAAGAAAAAGATCTTGGGACTTATAAAAAATCAACCCTTAAAACTGAAAAAATCACTCGCGGATTATTTTCAAACGATGAAATTACATTGATATATTTCTCAGAATATTCAAAAAGAATTGTCCAAGAGGTTTTTGTGTTTAATGTTGAAGACAAAAAAGTAAAATTAAAAGGTTATAGGTATGATTCGATTAATTAAATTTTTTGATGATGACTTCAATTAGTAGAACATGTATTTTTCGATAAAAGCACCGATGACTTTTTCATGAAGTTCTACAGAGCGCGAGAAACAGATGGTTCTGCGAGTCAGACGTTTGATTCAGGTTCTCAGGGGCAGATTCAGGGAAAGATCACAGGGAAAGATCATCGGTGTGCCGTAGGCACAAAGATGCAAAACGGGACGTTTTGCATGTAATAAGCATAAGCTTGTGAAGAGATGGTGGCAGGCCCACCGGGATCTGCTGTCGAAAGTCGGTCCCAAACAACCCCGGCTGCTGTGGTAAAGAGCCATGGTAGTGAGCGCAGGGGGAACGATCAGACGAGATCTGATCCTATGTGAGATTAGGGAAAGATCATCAAAAGTTCACTTTTTGTACTAAATAATTCGCATTTTATGTTTAAAAATTGAGATATTCCTCATTACCTGAAGCTGTTTTTTATTGCTTATACATGATCAAATACTCCTTACATAATTAAGGAGAACAAAATGGAACTTAAAAAATTGATGGAACATATTTCTATTATCCCCGATTACAGACAAGCCTGGAAAGTGGAACATAAATTATCGGATATTCTACTGTTGACTATTTGTGCCGTTATTTCTGGTGCAGAAGGCTGGGAAGATATAGAGGATTTTGGGGAAACACATCCCGATGTTTTGAAGTGATATGGTGATTTTGAAAATGGTATTCCTGTTCACGATATCATTGCCAGAGTTGTATCCTGTATCAGTCCTGCAAAATTTCACGAGTGCTTTATTAACTGGATGCGTGACTGCCATTCTTCAGATGATAAAGACGTCATTGCAATTGATGGAAAAACGCTCCGGCACTCTTATGACAAGAGTCGCCGCAGGGGAGCGATTCATGTCATTAGTGCGTTCTCAACAATGCACAGTCTGGTCATCGGACAGATCAAGACGGATGAGAAATCTAATGAGATTACAGCTATCCCAGAACTTCTTAACATGCTGGATATTAAAGGAAAAATCATCACAACTGATGCGATGGGTTGCCAGAAAGATATTGCAGAGAAGATACAAAAACAGGGAGGTGATTATTTATTCGCGGTAAAAGGAAACCAGGGGCGGCTAAATAAAGCCTTTGAGGAAAAATTTCCGCTGAAAGAATTAAATAATCCAGAGCATGACAGTTACGCAATTAGTGAAAAGAGTCACGGC",
    "contigLength": 5498578,
    "proteins": {
     "b0497": [
      "MSGKPAARQGDMTQYGGPIVQGSAGVRIGAPTGVACSVCPGGMTSGNPVNPLLGAKVLPGETDLALPGPLPFILSRTYSSYRTKTPAPVGVFGPGWKAPSDIRLQLRDDGLILNDNGGRSIHFEPLLPGEAVYSRSESMWLVRGGKAAQPDGHTLARLWGALPPDIRLSPHLYLATNSAQGPWWILGWSERVPGAEDVLPAPLPPYRVLTGMADRFGRTLTYRREAAGDLAGEITGVTDGAGREFRLVLTTQAQRAEEARTSSLSSSDSSRPLSASAFPDTLPGTEYGPDRGIRLSAVWLMHDPAYPESLPAAPLVRYTYTEAGELLAVYDRSNTQVRAFTYDAQHPGRMVAHRYAGRPEMRYRYDDTGRVVEQLNPAGLSYRYLYEQDRITVTDSLNRREVLHTEGGAGLKRVVKKELADGSVTRSGYDAAGRLTAQTDAAGRRTEYGLNVVSGDITDITTPDGRETKFYYNDGNQLTAVVSPDGLESRREYDEPGRLVSETSRSGETVRYRYDDAHSELPATTTDATGSTRQMTWSRYGQLLAFTDCSGYQTRYEYDRFGQMTAVHREEGISLYRRYDNRGRLTSVKDAQGRETRYEYNAAGDLTAVITPDGNRSETQYDAWGKAVSTTQGGLTRSMEYDAAGRVISLTNENGSHSVFSYDALDRLVQQGGFDGRTQRYHYDLTGKLTQSEDEGLVILWYYDESDRITHRTVNGEPAEQWQYDGHGWLTDISHLSEGHRVAVHYGYDDKGRLTGECQTVENPETGELLWQHETKHAYNEQGLANRVTPDSLPPVEWLTYGSGYLAGMKLGGTPLVEYTRDRLHRETVRSFGSMAGSNAAYELTSTYTPAGQLQSQHLNSLVYDRDYGWSDNGDLVRISGPRQTREYGYSATGRLESVRTLAPDLDIRIPYATDPAGNRLPDPELHPDSTLTVWPDNRIAEDAHYVYRHDEYGRLTEKTDRIPAGVIRTDDERTHHYHYDSQHRLVFYTRIQHGEPLVESRYLYDPLGRRMAKRVWRRERDLTGWMSLSRKPEVTWYGWDGDRLTTVQTDTTRIQTVYEPGSFTPLIRVETENGEREKAQRRSLAETLQQEGSENGHGVVFPAELVRLLDRLEEEIRADRVSSESRAWLAQCGLTVEQLARQVEPEYTPARKAHLYHCDHRGLPLALISEDGNTAWSAEYDEWGNQLNEENPHHVYQPYRLPGQQHDEESGLYYNRHRYYDPLQGRYITQDPMGLKGGWNLYQYPLNPLQQIDPMGLLQTWDDARSGACTGGVCGVLSRIIGPSKFDSTADAALDALKETQNRSLCNDMEYSGIVCKDTNGKYFASKAETDNLRKESYPLKRKCPTGTDRVAAYHTHGADSHGDYVDEFFSSSDKNLVRSKDNNLEAFYLATPDGRFEALNNKGEYIFIRNSVPGLSSVCIPYHD*",
      4281
     ]
    }
   }
  }
 }
};
