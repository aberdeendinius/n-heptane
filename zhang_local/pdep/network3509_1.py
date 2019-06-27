species(
    label = '[CH2]C(C=O)C(C)[O](12750)',
    structure = SMILES('[CH2]C(C=O)C(C)[O]'),
    E0 = (11.7569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,250.745,250.875],'cm^-1')),
        HinderedRotor(inertia=(0.00269962,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00269174,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25729,'amu*angstrom^2'), symmetry=1, barrier=(11.4093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254878,'amu*angstrom^2'), symmetry=1, barrier=(11.4127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907792,0.0720923,-7.6221e-05,4.68729e-08,-1.21665e-11,1521.9,27.3735], Tmin=(100,'K'), Tmax=(916.66,'K')), NASAPolynomial(coeffs=[9.34716,0.0352656,-1.59586e-05,3.04529e-09,-2.13453e-13,-25.3057,-12.6072], Tmin=(916.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.7569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC=O(606)',
    structure = SMILES('CC=O'),
    E0 = (-177.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1427.17,1427.17,1427.17,1427.17,3755.47],'cm^-1')),
        HinderedRotor(inertia=(0.717734,'amu*angstrom^2'), symmetry=1, barrier=(16.5021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70079,0.000387835,3.86929e-05,-4.52447e-08,1.58859e-11,-21380.9,9.13562], Tmin=(100,'K'), Tmax=(984.198,'K')), NASAPolynomial(coeffs=[4.58889,0.0128894,-4.91502e-06,9.26508e-10,-6.71011e-14,-22336,0.901072], Tmin=(984.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsH)"""),
)

species(
    label = 'C=CC=O(5269)',
    structure = SMILES('C=CC=O'),
    E0 = (-81.3387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.873408,'amu*angstrom^2'), symmetry=1, barrier=(20.0814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3136.31,'J/mol'), sigma=(5.14154,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.88 K, Pc=52.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.9738,0.0193269,-1.02836e-06,-7.40922e-09,2.6466e-12,-9743.32,12.1361], Tmin=(100,'K'), Tmax=(1315.19,'K')), NASAPolynomial(coeffs=[7.40832,0.0154746,-7.62321e-06,1.50372e-09,-1.06406e-13,-11743,-13.6408], Tmin=(1315.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.3387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC([O])C1CC1[O](13681)',
    structure = SMILES('CC([O])C1CC1[O]'),
    E0 = (101.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14857,0.0469664,1.22132e-05,-5.27288e-08,2.43372e-11,12294.6,24.9957], Tmin=(100,'K'), Tmax=(968.687,'K')), NASAPolynomial(coeffs=[16.5242,0.0215622,-7.42533e-06,1.37547e-09,-1.01306e-13,7528.78,-57.9174], Tmin=(968.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1C(C)OC1[O](13664)',
    structure = SMILES('[CH2]C1C(C)OC1[O]'),
    E0 = (64.9921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66912,0.0426014,9.69822e-06,-4.43613e-08,2.21015e-11,7908.9,24.1228], Tmin=(100,'K'), Tmax=(876.597,'K')), NASAPolynomial(coeffs=[10.3078,0.027971,-7.68453e-06,1.11797e-09,-6.90892e-14,5441.96,-21.8486], Tmin=(876.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.9921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = 'H(8)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(C=O)C(C)[O](13682)',
    structure = SMILES('C=C(C=O)C(C)[O]'),
    E0 = (-89.7735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,452.943,465.3],'cm^-1')),
        HinderedRotor(inertia=(0.0872359,'amu*angstrom^2'), symmetry=1, barrier=(13.1277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.083846,'amu*angstrom^2'), symmetry=1, barrier=(13.1244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0865967,'amu*angstrom^2'), symmetry=1, barrier=(13.124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12556,0.0580179,-4.1184e-05,1.38435e-08,-1.85252e-12,-10689.8,26.111], Tmin=(100,'K'), Tmax=(1729.11,'K')), NASAPolynomial(coeffs=[16.4797,0.0224982,-1.03702e-05,1.96299e-09,-1.3477e-13,-15999.6,-56.3717], Tmin=(1729.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.7735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(C=O)C(C)=O(12664)',
    structure = SMILES('[CH2]C(C=O)C(C)=O'),
    E0 = (-159.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,1128.59],'cm^-1')),
        HinderedRotor(inertia=(0.00645083,'amu*angstrom^2'), symmetry=1, barrier=(5.87153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161744,'amu*angstrom^2'), symmetry=1, barrier=(5.88879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169154,'amu*angstrom^2'), symmetry=1, barrier=(5.91846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493225,'amu*angstrom^2'), symmetry=1, barrier=(18.3543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56472,0.0546489,-3.96865e-05,1.51477e-08,-2.45314e-12,-19123.9,25.8247], Tmin=(100,'K'), Tmax=(1382.74,'K')), NASAPolynomial(coeffs=[9.88305,0.0305858,-1.35831e-05,2.56243e-09,-1.77746e-13,-21424.4,-17.0023], Tmin=(1382.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C[CH][O](605)',
    structure = SMILES('C[CH][O]'),
    E0 = (149.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2066.51],'cm^-1')),
        HinderedRotor(inertia=(0.362113,'amu*angstrom^2'), symmetry=1, barrier=(8.32568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.20363,0.021847,-3.14755e-05,3.43227e-08,-1.42322e-11,17997,11.0861], Tmin=(100,'K'), Tmax=(846.374,'K')), NASAPolynomial(coeffs=[1.2024,0.020386,-9.53523e-06,1.79858e-09,-1.23081e-13,18726.8,22.7175], Tmin=(846.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH]=O(373)',
    structure = SMILES('[CH]=O'),
    E0 = (33.3613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([880.033,2045.55,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.35602,-0.00347087,1.25664e-05,-9.99488e-09,2.27888e-12,3995.77,2.75113], Tmin=(100,'K'), Tmax=(1565.72,'K')), NASAPolynomial(coeffs=[4.61867,0.00504456,-4.39241e-06,9.73282e-10,-7.07436e-14,2787.5,-2.22961], Tmin=(1565.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH) + radical(HCdsJO)"""),
)

species(
    label = 'C=CC(C)[O](738)',
    structure = SMILES('C=CC(C)[O]'),
    E0 = (41.0323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,460.346,460.818],'cm^-1')),
        HinderedRotor(inertia=(0.0878022,'amu*angstrom^2'), symmetry=1, barrier=(13.2705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883523,'amu*angstrom^2'), symmetry=1, barrier=(13.2736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10805,0.0309456,1.3595e-05,-3.86994e-08,1.66088e-11,5012.63,19.9443], Tmin=(100,'K'), Tmax=(988.834,'K')), NASAPolynomial(coeffs=[11.3016,0.0192069,-7.20555e-06,1.35312e-09,-9.81465e-14,1950.16,-30.5978], Tmin=(988.834,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.0323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=C[O](5266)',
    structure = SMILES('[CH2]C=C[O]'),
    E0 = (90.2929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.57685,'amu*angstrom^2'), symmetry=1, barrier=(36.2549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.69019,0.0144913,4.15491e-05,-7.27602e-08,3.14101e-11,10920.2,13.4175], Tmin=(100,'K'), Tmax=(922.751,'K')), NASAPolynomial(coeffs=[14.044,0.00224417,1.35973e-06,-3.04875e-10,1.62832e-14,7250.86,-48.974], Tmin=(922.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.2929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH3](11)',
    structure = SMILES('[CH3]'),
    E0 = (135.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([570.572,1408.13,1408.49,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91547,0.00184154,3.48742e-06,-3.32748e-09,8.49957e-13,16285.6,0.351741], Tmin=(100,'K'), Tmax=(1337.63,'K')), NASAPolynomial(coeffs=[3.54146,0.00476787,-1.82148e-06,3.28877e-10,-2.22546e-14,16224,1.66035], Tmin=(1337.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""),
)

species(
    label = '[CH2]C(C=O)C=O(12779)',
    structure = SMILES('[CH2]C(C=O)C=O'),
    E0 = (-104.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.36304,'amu*angstrom^2'), symmetry=1, barrier=(31.339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000540802,'amu*angstrom^2'), symmetry=1, barrier=(4.85933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36291,'amu*angstrom^2'), symmetry=1, barrier=(31.336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09608,0.0407753,-2.78866e-05,9.59164e-09,-1.36393e-12,-12541.8,21.6355], Tmin=(100,'K'), Tmax=(1580.53,'K')), NASAPolynomial(coeffs=[9.87979,0.0210762,-9.19104e-06,1.7058e-09,-1.16576e-13,-15002.3,-19.4794], Tmin=(1580.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CC(=C[O])C(C)[O](13683)',
    structure = SMILES('CC(=C[O])C(C)[O]'),
    E0 = (-65.3528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,341.619,341.619,341.619],'cm^-1')),
        HinderedRotor(inertia=(0.202267,'amu*angstrom^2'), symmetry=1, barrier=(16.7508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202267,'amu*angstrom^2'), symmetry=1, barrier=(16.7508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202266,'amu*angstrom^2'), symmetry=1, barrier=(16.7508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.645401,0.0596775,-1.79466e-05,-2.65285e-08,1.63198e-11,-7726.63,27.2194], Tmin=(100,'K'), Tmax=(966.321,'K')), NASAPolynomial(coeffs=[18.5905,0.0187667,-6.24314e-06,1.13529e-09,-8.31109e-14,-12752.8,-66.8023], Tmin=(966.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.3528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(C=O)[C](C)O(13684)',
    structure = SMILES('[CH2]C(C=O)[C](C)O'),
    E0 = (-41.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.332726,0.0880793,-0.000131657,1.12039e-07,-3.78706e-11,-4923.62,29.0395], Tmin=(100,'K'), Tmax=(840.337,'K')), NASAPolynomial(coeffs=[8.64145,0.0360691,-1.6577e-05,3.09678e-09,-2.10707e-13,-5880.08,-6.98215], Tmin=(840.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C[C]([O])C(C)C=O(13685)',
    structure = SMILES('C[C]([O])C(C)C=O'),
    E0 = (-22.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949662,0.0718974,-8.1487e-05,5.68673e-08,-1.69697e-11,-2555.64,26.2951], Tmin=(100,'K'), Tmax=(800.476,'K')), NASAPolynomial(coeffs=[7.82273,0.0375515,-1.71251e-05,3.26288e-09,-2.27819e-13,-3655.96,-5.33366], Tmin=(800.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])C(C)[C]=O(13686)',
    structure = SMILES('CC([O])C(C)[C]=O'),
    E0 = (-40.0652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08166,0.0629516,-5.09272e-05,2.21627e-08,-4.01898e-12,-4712.79,27.5396], Tmin=(100,'K'), Tmax=(1282.62,'K')), NASAPolynomial(coeffs=[11.5128,0.0304208,-1.2883e-05,2.38842e-09,-1.64689e-13,-7388.62,-25.3809], Tmin=(1282.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.0652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C(=C[O])C(C)O(13687)',
    structure = SMILES('[CH2]C(=C[O])C(C)O'),
    E0 = (-144.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499589,0.0609652,-1.45013e-05,-3.60057e-08,2.12808e-11,-17204.3,28.2155], Tmin=(100,'K'), Tmax=(946.941,'K')), NASAPolynomial(coeffs=[20.8059,0.0147468,-3.95132e-06,6.82384e-10,-5.198e-14,-22823.7,-78.0087], Tmin=(946.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)C([CH2])C=O(12792)',
    structure = SMILES('[CH2]C(O)C([CH2])C=O'),
    E0 = (-7.01478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.428534,0.0833066,-0.000109453,8.16414e-08,-2.48125e-11,-719.471,29.9127], Tmin=(100,'K'), Tmax=(800.968,'K')), NASAPolynomial(coeffs=[10.5636,0.0326932,-1.46686e-05,2.75145e-09,-1.89506e-13,-2343.07,-16.734], Tmin=(800.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.01478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([C]=O)C(C)O(13688)',
    structure = SMILES('[CH2]C([C]=O)C(C)O'),
    E0 = (-59.9152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728736,0.0756568,-8.73348e-05,5.73964e-08,-1.55558e-11,-7091.57,29.3623], Tmin=(100,'K'), Tmax=(889.8,'K')), NASAPolynomial(coeffs=[10.3606,0.0323557,-1.43355e-05,2.70037e-09,-1.87583e-13,-8805.58,-15.9809], Tmin=(889.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.9152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)CJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([O])C(C)C=O(13689)',
    structure = SMILES('[CH2]C([O])C(C)C=O'),
    E0 = (12.8353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0133696,'amu*angstrom^2'), symmetry=1, barrier=(2.10586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133493,'amu*angstrom^2'), symmetry=1, barrier=(2.06266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761204,'amu*angstrom^2'), symmetry=1, barrier=(11.6816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507167,'amu*angstrom^2'), symmetry=1, barrier=(11.6608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.996761,0.0679205,-6.30765e-05,3.27318e-08,-7.13073e-12,1650.32,27.3283], Tmin=(100,'K'), Tmax=(1080.86,'K')), NASAPolynomial(coeffs=[10.4702,0.0328623,-1.44243e-05,2.72406e-09,-1.90165e-13,-397.611,-19.1124], Tmin=(1080.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.8353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][CH]C(C)[O](730)',
    structure = SMILES('[CH2][CH]C(C)[O]'),
    E0 = (318.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,717.471,717.473],'cm^-1')),
        HinderedRotor(inertia=(0.10698,'amu*angstrom^2'), symmetry=1, barrier=(2.45968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10698,'amu*angstrom^2'), symmetry=1, barrier=(2.45968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00673363,'amu*angstrom^2'), symmetry=1, barrier=(2.45968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96464,0.0376728,-1.08566e-05,-1.06109e-08,6.18527e-12,38414.4,23.3594], Tmin=(100,'K'), Tmax=(1031.8,'K')), NASAPolynomial(coeffs=[10.0599,0.0210439,-8.13097e-06,1.48667e-09,-1.03787e-13,35958.5,-19.7546], Tmin=(1031.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(=C[O])C(C)[O](13690)',
    structure = SMILES('[CH2]C(=C[O])C(C)[O]'),
    E0 = (86.1464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,400.33,400.344,400.351],'cm^-1')),
        HinderedRotor(inertia=(0.197563,'amu*angstrom^2'), symmetry=1, barrier=(22.4685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197571,'amu*angstrom^2'), symmetry=1, barrier=(22.4686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197564,'amu*angstrom^2'), symmetry=1, barrier=(22.4687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659317,0.0570082,-8.72707e-06,-4.08642e-08,2.28079e-11,10496.3,27.4503], Tmin=(100,'K'), Tmax=(948.261,'K')), NASAPolynomial(coeffs=[20.7476,0.0129429,-3.35933e-06,5.93486e-10,-4.68803e-14,4857.93,-78.0388], Tmin=(948.261,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.1464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH][O])C=O(13547)',
    structure = SMILES('[CH2]C([CH][O])C=O'),
    E0 = (223.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,1656.02],'cm^-1')),
        HinderedRotor(inertia=(0.236434,'amu*angstrom^2'), symmetry=1, barrier=(5.43608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237351,'amu*angstrom^2'), symmetry=1, barrier=(5.45718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238528,'amu*angstrom^2'), symmetry=1, barrier=(5.48423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2747,0.0721842,-0.000133189,1.30951e-07,-4.80698e-11,27006,24.7433], Tmin=(100,'K'), Tmax=(867.505,'K')), NASAPolynomial(coeffs=[3.3884,0.0336322,-1.67212e-05,3.1697e-09,-2.1486e-13,27723.2,21.0936], Tmin=(867.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(C=O)[C](C)[O](13691)',
    structure = SMILES('[CH2]C(C=O)[C](C)[O]'),
    E0 = (188.385,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,360,370,350,254.889,254.908],'cm^-1')),
        HinderedRotor(inertia=(0.00259441,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150281,'amu*angstrom^2'), symmetry=1, barrier=(6.92928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150296,'amu*angstrom^2'), symmetry=1, barrier=(6.92921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150301,'amu*angstrom^2'), symmetry=1, barrier=(6.92916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497708,0.0840636,-0.000125704,1.07012e-07,-3.63143e-11,22776.8,28.2554], Tmin=(100,'K'), Tmax=(830.485,'K')), NASAPolynomial(coeffs=[8.53482,0.0343487,-1.60338e-05,3.0195e-09,-2.06578e-13,21821.4,-6.74141], Tmin=(830.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([O])C([CH2])C=O(12644)',
    structure = SMILES('[CH2]C([O])C([CH2])C=O'),
    E0 = (223.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,237.377,2887.88],'cm^-1')),
        HinderedRotor(inertia=(0.31931,'amu*angstrom^2'), symmetry=1, barrier=(12.7646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00299155,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319235,'amu*angstrom^2'), symmetry=1, barrier=(12.7648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788709,'amu*angstrom^2'), symmetry=1, barrier=(31.5423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.624896,0.0788821,-0.000101868,7.41622e-08,-2.20348e-11,26979.6,29.0184], Tmin=(100,'K'), Tmax=(817.816,'K')), NASAPolynomial(coeffs=[10.399,0.0310773,-1.41882e-05,2.68947e-09,-1.86674e-13,25380.9,-16.1705], Tmin=(817.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CC(C)OJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([C]=O)C(C)[O](13692)',
    structure = SMILES('[CH2]C([C]=O)C(C)[O]'),
    E0 = (170.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,465.63,465.631],'cm^-1')),
        HinderedRotor(inertia=(0.0600935,'amu*angstrom^2'), symmetry=1, barrier=(9.24578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0600945,'amu*angstrom^2'), symmetry=1, barrier=(9.24581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0600934,'amu*angstrom^2'), symmetry=1, barrier=(9.24578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0600933,'amu*angstrom^2'), symmetry=1, barrier=(9.2458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.914242,0.0713934,-8.04757e-05,5.11034e-08,-1.33997e-11,20608,28.5048], Tmin=(100,'K'), Tmax=(917.113,'K')), NASAPolynomial(coeffs=[10.2639,0.0306144,-1.37784e-05,2.6195e-09,-1.83134e-13,18893,-15.7929], Tmin=(917.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CC(C)OJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CC([O])C1[CH]OC1(13693)',
    structure = SMILES('CC([O])C1[CH]OC1'),
    E0 = (85.4249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0942,0.0515592,2.77384e-06,-5.17409e-08,2.85818e-11,10390.7,23.7291], Tmin=(100,'K'), Tmax=(871.295,'K')), NASAPolynomial(coeffs=[16.2049,0.0184363,-2.60735e-06,1.25026e-10,-4.76473e-16,6381.6,-54.9855], Tmin=(871.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.4249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCsJOCs) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1[CH]OOC1C(13618)',
    structure = SMILES('[CH2]C1[CH]OOC1C'),
    E0 = (197.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42202,0.041745,3.09707e-05,-8.06206e-08,3.89465e-11,23903.9,22.7662], Tmin=(100,'K'), Tmax=(873.311,'K')), NASAPolynomial(coeffs=[15.6824,0.0186917,-2.02404e-06,-1.87941e-11,9.55312e-15,19801.5,-53.3272], Tmin=(873.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(Isobutyl) + radical(CCsJOOC)"""),
)

species(
    label = 'C=C(C=O)C(C)O(13694)',
    structure = SMILES('C=C(C=O)C(C)O'),
    E0 = (-320.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06664,0.0609153,-4.37886e-05,1.52398e-08,-2.13785e-12,-38395.2,26.505], Tmin=(100,'K'), Tmax=(1638.11,'K')), NASAPolynomial(coeffs=[15.3432,0.0260541,-1.18664e-05,2.24819e-09,-1.55123e-13,-43072.5,-49.4171], Tmin=(1638.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'CC(=O)C(C)C=O(13695)',
    structure = SMILES('CC(=O)C(C)C=O'),
    E0 = (-370.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36979,0.0500484,-2.18055e-05,-7.11311e-10,2.05988e-12,-44427.9,26.1914], Tmin=(100,'K'), Tmax=(1238.93,'K')), NASAPolynomial(coeffs=[11.2399,0.0303742,-1.2747e-05,2.35761e-09,-1.62248e-13,-47809.3,-27.3171], Tmin=(1238.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-370.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH)"""),
)

species(
    label = '[CH2][C](C[O])C(C)[O](9641)',
    structure = SMILES('[CH2][C](C[O])C(C)[O]'),
    E0 = (307.646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,839.257,840.213,842.082],'cm^-1')),
        HinderedRotor(inertia=(0.178196,'amu*angstrom^2'), symmetry=1, barrier=(4.09707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179396,'amu*angstrom^2'), symmetry=1, barrier=(4.12467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00818666,'amu*angstrom^2'), symmetry=1, barrier=(4.1081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0083241,'amu*angstrom^2'), symmetry=1, barrier=(4.15026,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31278,0.0628814,-5.83521e-05,3.52611e-08,-9.59684e-12,37094.7,28.9349], Tmin=(100,'K'), Tmax=(853.753,'K')), NASAPolynomial(coeffs=[6.29143,0.0395556,-1.73702e-05,3.2599e-09,-2.26155e-13,36244.6,5.70296], Tmin=(853.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(C[O])[C](C)[O](13696)',
    structure = SMILES('[CH2]C(C[O])[C](C)[O]'),
    E0 = (331.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,244.271,1912.18,1912.49,1912.81],'cm^-1')),
        HinderedRotor(inertia=(1.48989,'amu*angstrom^2'), symmetry=1, barrier=(63.0211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177076,'amu*angstrom^2'), symmetry=1, barrier=(7.54594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176636,'amu*angstrom^2'), symmetry=1, barrier=(7.54419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177611,'amu*angstrom^2'), symmetry=1, barrier=(7.54426,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676731,0.0809446,-0.000117113,1.02178e-07,-3.54321e-11,40006.4,30.4337], Tmin=(100,'K'), Tmax=(855.681,'K')), NASAPolynomial(coeffs=[5.98888,0.0402616,-1.80106e-05,3.31848e-09,-2.23864e-13,39677.5,9.02431], Tmin=(855.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]O)[C](C)[O](13697)',
    structure = SMILES('[CH2]C([CH]O)[C](C)[O]'),
    E0 = (286.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0638016,0.0922029,-0.000135157,1.0808e-07,-3.38422e-11,34569.4,31.3225], Tmin=(100,'K'), Tmax=(890.132,'K')), NASAPolynomial(coeffs=[11.2542,0.0307225,-1.26901e-05,2.23085e-09,-1.45677e-13,33020.7,-18.8713], Tmin=(890.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCsJOH) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C([CH2])C[O](13698)',
    structure = SMILES('[CH2]C([O])C([CH2])C[O]'),
    E0 = (366.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,990.419,990.451,990.48],'cm^-1')),
        HinderedRotor(inertia=(0.155487,'amu*angstrom^2'), symmetry=1, barrier=(3.57496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155361,'amu*angstrom^2'), symmetry=1, barrier=(3.57206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15546,'amu*angstrom^2'), symmetry=1, barrier=(3.57434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155513,'amu*angstrom^2'), symmetry=1, barrier=(3.57555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712346,0.0769615,-9.80807e-05,7.65611e-08,-2.47525e-11,44213.1,31.5178], Tmin=(100,'K'), Tmax=(816.54,'K')), NASAPolynomial(coeffs=[8.05093,0.0366322,-1.59493e-05,2.93585e-09,-1.99492e-13,43160.6,-1.50515], Tmin=(816.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([O])C([CH2])[CH]O(13699)',
    structure = SMILES('[CH2]C([O])C([CH2])[CH]O'),
    E0 = (321.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.22839,0.0864709,-0.000108747,7.06188e-08,-1.68173e-11,38770.6,31.9575], Tmin=(100,'K'), Tmax=(743.957,'K')), NASAPolynomial(coeffs=[13.1897,0.0273313,-1.07761e-05,1.88477e-09,-1.2445e-13,36550.2,-28.7015], Tmin=(743.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCsJOH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'CH2(S)(14)',
    structure = SMILES('[CH2]'),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144068,5.45069e-06,-3.58002e-09,7.56192e-13,50400.6,-0.411765], Tmin=(100,'K'), Tmax=(1442.36,'K')), NASAPolynomial(coeffs=[2.62648,0.00394763,-1.49924e-06,2.54539e-10,-1.62956e-14,50691.8,6.78378], Tmin=(1442.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C(C=O)C[O](12622)',
    structure = SMILES('[CH2]C(C=O)C[O]'),
    E0 = (43.5244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,1904.17],'cm^-1')),
        HinderedRotor(inertia=(0.220423,'amu*angstrom^2'), symmetry=1, barrier=(5.06795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220438,'amu*angstrom^2'), symmetry=1, barrier=(5.06831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220384,'amu*angstrom^2'), symmetry=1, barrier=(5.06706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3870.7,'J/mol'), sigma=(6.40089,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.59 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52783,0.063651,-0.000102817,9.9992e-08,-3.76092e-11,5314.83,24.2849], Tmin=(100,'K'), Tmax=(842.416,'K')), NASAPolynomial(coeffs=[2.79935,0.0367788,-1.78705e-05,3.40907e-09,-2.34164e-13,5839.88,22.7565], Tmin=(842.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.5244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[C-]#[O+](374)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (299.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33667,0.00896487,-2.66756e-05,3.61071e-08,-1.57199e-11,36069.2,-1.20266], Tmin=(100,'K'), Tmax=(865.594,'K')), NASAPolynomial(coeffs=[-0.394107,0.0117562,-6.47408e-06,1.26375e-09,-8.67562e-14,37256.3,19.3844], Tmin=(865.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.89,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH2]CC(C)[O](3574)',
    structure = SMILES('[CH2]CC(C)[O]'),
    E0 = (118.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,792.118],'cm^-1')),
        HinderedRotor(inertia=(0.516548,'amu*angstrom^2'), symmetry=1, barrier=(11.8765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125741,'amu*angstrom^2'), symmetry=1, barrier=(2.89104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0063771,'amu*angstrom^2'), symmetry=1, barrier=(2.85618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3512.21,'J/mol'), sigma=(6.17291,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.60 K, Pc=33.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7348,0.0432476,-2.04579e-05,-3.21265e-10,2.22143e-12,14379.3,21.6911], Tmin=(100,'K'), Tmax=(1125.04,'K')), NASAPolynomial(coeffs=[10.0238,0.0244802,-9.70653e-06,1.76441e-09,-1.21234e-13,11836.8,-22.2855], Tmin=(1125.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])CC=C[O](12748)',
    structure = SMILES('CC([O])CC=C[O]'),
    E0 = (-49.7417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.717566,0.0585811,-1.68609e-05,-2.62176e-08,1.59217e-11,-5852.1,27.5392], Tmin=(100,'K'), Tmax=(966.356,'K')), NASAPolynomial(coeffs=[17.8565,0.0198209,-6.65052e-06,1.20073e-09,-8.69923e-14,-10667.2,-62.3347], Tmin=(966.356,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.7417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC([O])[CH]CC=O(13700)',
    structure = SMILES('CC([O])[CH]CC=O'),
    E0 = (6.90878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35907,0.0593895,-4.45329e-05,1.80466e-08,-3.12792e-12,924.654,28.0197], Tmin=(100,'K'), Tmax=(1301.03,'K')), NASAPolynomial(coeffs=[9.7709,0.0335276,-1.47161e-05,2.76811e-09,-1.92089e-13,-1264.17,-14.7762], Tmin=(1301.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.90878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'CC1OCC1C=O(12758)',
    structure = SMILES('CC1OCC1C=O'),
    E0 = (-245.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64292,0.0415241,1.49025e-05,-5.09248e-08,2.44592e-11,-29382.1,21.2803], Tmin=(100,'K'), Tmax=(891.681,'K')), NASAPolynomial(coeffs=[11.313,0.0263772,-7.10963e-06,1.04047e-09,-6.56144e-14,-32229,-30.5572], Tmin=(891.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-245.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]C=COC(C)[O](12749)',
    structure = SMILES('[CH2]C=COC(C)[O]'),
    E0 = (-44.5593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,306.832,306.837,306.84,306.841],'cm^-1')),
        HinderedRotor(inertia=(0.359504,'amu*angstrom^2'), symmetry=1, barrier=(24.0189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359506,'amu*angstrom^2'), symmetry=1, barrier=(24.0189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359491,'amu*angstrom^2'), symmetry=1, barrier=(24.0188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00179046,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544477,0.063415,-3.02216e-05,-1.1181e-08,9.96948e-12,-5223.7,27.194], Tmin=(100,'K'), Tmax=(997.596,'K')), NASAPolynomial(coeffs=[17.9064,0.0212082,-7.97061e-06,1.48994e-09,-1.0766e-13,-10051.6,-63.3613], Tmin=(997.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.5593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=CO)C(C)[O](13701)',
    structure = SMILES('[CH2]C(=CO)C(C)[O]'),
    E0 = (-55.3162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,219.929,219.943],'cm^-1')),
        HinderedRotor(inertia=(0.616441,'amu*angstrom^2'), symmetry=1, barrier=(21.1588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616388,'amu*angstrom^2'), symmetry=1, barrier=(21.158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616376,'amu*angstrom^2'), symmetry=1, barrier=(21.1589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616314,'amu*angstrom^2'), symmetry=1, barrier=(21.1589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246548,0.0639887,-1.21373e-05,-4.71322e-08,2.76099e-11,-6500.7,27.522], Tmin=(100,'K'), Tmax=(928.227,'K')), NASAPolynomial(coeffs=[24.0311,0.00920481,-7.06617e-07,3.1833e-11,-6.71272e-15,-12971.6,-96.5245], Tmin=(928.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.3162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CC(C)OJ)"""),
)

species(
    label = 'CH2(T)(28)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.58,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CC([O])C=C[O](1111)',
    structure = SMILES('CC([O])C=C[O]'),
    E0 = (-26.2977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,273.512,274.01,274.198],'cm^-1')),
        HinderedRotor(inertia=(0.38108,'amu*angstrom^2'), symmetry=1, barrier=(20.3493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382682,'amu*angstrom^2'), symmetry=1, barrier=(20.3508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48963,0.0396646,1.26669e-05,-5.3726e-08,2.58057e-11,-3058.29,23.4328], Tmin=(100,'K'), Tmax=(945.207,'K')), NASAPolynomial(coeffs=[17.6459,0.0102077,-2.34143e-06,4.16246e-10,-3.49302e-14,-7850.82,-62.7968], Tmin=(945.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.2977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O(T)(63)',
    structure = SMILES('[O]'),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29230.2,4.09104], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29230.2,4.09104], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C([CH]C)C=O(12927)',
    structure = SMILES('[CH2]C([CH]C)C=O'),
    E0 = (153.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,298.079],'cm^-1')),
        HinderedRotor(inertia=(0.0937562,'amu*angstrom^2'), symmetry=1, barrier=(5.91007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00189972,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00190097,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0937914,'amu*angstrom^2'), symmetry=1, barrier=(5.90927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62238,0.0540135,-4.22432e-05,1.8639e-08,-3.56672e-12,18578.2,25.1257], Tmin=(100,'K'), Tmax=(1186.45,'K')), NASAPolynomial(coeffs=[8.33854,0.0313704,-1.36159e-05,2.55323e-09,-1.77214e-13,16984.5,-8.42406], Tmin=(1186.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]C(C=O)C(C)[O](13702)',
    structure = SMILES('[CH]C(C=O)C(C)[O]'),
    E0 = (249.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09515,0.0657977,-6.28544e-05,3.32363e-08,-7.33645e-12,30106.3,27.3459], Tmin=(100,'K'), Tmax=(1070.93,'K')), NASAPolynomial(coeffs=[10.5064,0.0306453,-1.36171e-05,2.58489e-09,-1.80968e-13,28090.6,-18.7027], Tmin=(1070.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (11.7569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (101.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (133.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (132.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (123.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (77.3144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (94.0545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (11.7569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (62.5571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (137.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (177.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (129.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (129.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (87.0376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (150.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (139.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (95.6785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (239.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (352.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (297.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (359.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (400.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (435.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (382.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (198.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (197.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (75.1571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (75.1571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (329.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (395.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (294.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (391.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (346.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (462.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (724.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (171.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (134.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (20.0412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (269.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (88.5506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (389.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (560.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (461.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['CC=O(606)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['CC([O])C1CC1[O](13681)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.93521e+09,'s^-1'), n=0.743095, Ea=(89.4945,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 85.9 to 89.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['[CH2]C1C(C)OC1[O](13664)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C=C(C=O)C(C)[O](13682)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(72.1434,'m^3/(mol*s)'), n=1.66666, Ea=(10.8177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;HJ] for rate rule [Cds-COCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH2]C(C=O)C(C)=O(12664)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C[CH][O](605)', 'C=CC=O(5269)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.76206,'m^3/(mol*s)'), n=1.97634, Ea=(9.22116,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-OneDeH_Cds;CJ] + [Cds-COH_Cds;YJ] for rate rule [Cds-COH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=O(373)', 'C=CC(C)[O](738)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CJ] for rate rule [Cds-CsH_Cds-HH;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC=O(606)', '[CH2]C=C[O](5266)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(99.3703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 98.2 to 99.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH3](11)', '[CH2]C(C=O)C=O(12779)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.22517,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-CsH_O;CsJ-HHH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC(=C[O])C(C)[O](13683)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['[CH2]C(C=O)[C](C)O(13684)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['C[C]([O])C(C)C=O(13685)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['CC([O])C(C)[C]=O(13686)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_2H;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['[CH2]C(=C[O])C(C)O(13687)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['[CH2]C([C]=O)C(C)O(13688)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.75172e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([O])C(C)C=O(13689)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH][O](605)', '[CH2]C=C[O](5266)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=O(373)', '[CH2][CH]C(C)[O](730)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.7839e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[CH2]C(=C[O])C(C)[O](13690)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH3](11)', '[CH2]C([CH][O])C=O(13547)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH2]C(C=O)[C](C)[O](13691)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2]C([O])C([CH2])C=O(12644)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH2]C([C]=O)C(C)[O](13692)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['CC([O])C1[CH]OC1(13693)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['[CH2]C1[CH]OOC1C(13618)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(186.101,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 182.1 to 186.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['C=C(C=O)C(C)O(13694)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['CC(=O)C(C)C=O(13695)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C](C[O])C(C)[O](9641)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C[O])[C](C)[O](13696)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([CH]O)[C](C)[O](13697)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([O])C([CH2])C[O](13698)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])C([CH2])[CH]O(13699)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CH2(S)(14)', '[CH2]C(C=O)C[O](12622)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['[C-]#[O+](374)', '[CH2]CC(C)[O](3574)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['CC([O])CC=C[O](12748)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['CC([O])[CH]CC=O(13700)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['CC1OCC1C=O(12758)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C=COC(C)[O](12749)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(=CO)C(C)[O](13701)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(T)(28)', 'CC([O])C=C[O](1111)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O(T)(63)', '[CH2]C([CH]C)C=O(12927)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(8)', '[CH]C(C=O)C(C)[O](13702)'],
    products = ['[CH2]C(C=O)C(C)[O](12750)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '3509',
    isomers = [
        '[CH2]C(C=O)C(C)[O](12750)',
    ],
    reactants = [
        ('CC=O(606)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3509',
    Tmin = (1200,'K'),
    Tmax = (1500,'K'),
    Tcount = 10,
    Tlist = ([1201.48,1213.22,1236.21,1269.31,1310.55,1356.92,1404.16,1447.02,1479.84,1497.7],'K'),
    Pmin = (1,'atm'),
    Pmax = (10,'atm'),
    Pcount = 10,
    Plist = ([1.02771,1.14872,1.41959,1.89986,2.67608,3.83649,5.40396,7.23219,8.93758,9.98989],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

