species(
    label = 'C=[C]C(O)C(C)[O](13753)',
    structure = SMILES('C=[C]C(O)C(C)[O]'),
    E0 = (82.6234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,386.026,386.617,388.566],'cm^-1')),
        HinderedRotor(inertia=(0.00111578,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110548,'amu*angstrom^2'), symmetry=1, barrier=(11.7223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111143,'amu*angstrom^2'), symmetry=1, barrier=(11.7451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110928,'amu*angstrom^2'), symmetry=1, barrier=(11.7536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.697401,0.0679222,-5.9482e-05,2.71022e-08,-4.97994e-12,10060.1,30.0851], Tmin=(100,'K'), Tmax=(1299.39,'K')), NASAPolynomial(coeffs=[14.5232,0.0253612,-1.03503e-05,1.89464e-09,-1.30057e-13,6467.03,-40.2374], Tmin=(1299.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.6234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S)"""),
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
    label = 'C=C=CO(12571)',
    structure = SMILES('C=C=CO'),
    E0 = (-26.0646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.34368,'amu*angstrom^2'), symmetry=1, barrier=(30.8938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3437.21,'J/mol'), sigma=(5.57865,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.88 K, Pc=44.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31583,0.0236137,2.05754e-05,-5.73733e-08,2.79863e-11,-3061.58,12.125], Tmin=(100,'K'), Tmax=(901.949,'K')), NASAPolynomial(coeffs=[16.2977,-0.00239911,3.975e-06,-8.57293e-10,5.72973e-14,-7047.88,-62.0029], Tmin=(901.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.0646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C=[C]C(O)C(C)=O(27722)',
    structure = SMILES('C=[C]C(O)C(C)=O'),
    E0 = (-87.3665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,1685,370,1380,1390,370,380,2900,435,335.61,335.626],'cm^-1')),
        HinderedRotor(inertia=(0.00149658,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00149706,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116063,'amu*angstrom^2'), symmetry=1, barrier=(9.27398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116048,'amu*angstrom^2'), symmetry=1, barrier=(9.27393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26888,0.0574432,-4.53322e-05,1.85772e-08,-3.1259e-12,-10407.5,28.1243], Tmin=(100,'K'), Tmax=(1384.49,'K')), NASAPolynomial(coeffs=[12.2349,0.0257608,-1.10067e-05,2.04874e-09,-1.41343e-13,-13444,-28.3484], Tmin=(1384.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.3665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=C=C(O)C(C)[O](28168)',
    structure = SMILES('C=C=C(O)C(C)[O]'),
    E0 = (-32.9476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.797876,'amu*angstrom^2'), symmetry=1, barrier=(18.3447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799273,'amu*angstrom^2'), symmetry=1, barrier=(18.3769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.798423,'amu*angstrom^2'), symmetry=1, barrier=(18.3573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491931,0.0650016,-3.90332e-05,-6.33726e-09,9.95401e-12,-3825.36,26.6612], Tmin=(100,'K'), Tmax=(958.209,'K')), NASAPolynomial(coeffs=[19.5686,0.0140817,-4.27289e-06,7.52602e-10,-5.52773e-14,-8799.5,-71.4371], Tmin=(958.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.9476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC(O)C(C)[O](28169)',
    structure = SMILES('C#CC(O)C(C)[O]'),
    E0 = (10.7317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,483.613,483.921],'cm^-1')),
        HinderedRotor(inertia=(0.212365,'amu*angstrom^2'), symmetry=1, barrier=(35.2644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.606131,'amu*angstrom^2'), symmetry=1, barrier=(13.9361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184579,'amu*angstrom^2'), symmetry=1, barrier=(13.9359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.083941,'amu*angstrom^2'), symmetry=1, barrier=(13.9366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.345565,0.0704448,-6.85821e-05,3.42922e-08,-6.68929e-12,1430.72,28.4914], Tmin=(100,'K'), Tmax=(1296.78,'K')), NASAPolynomial(coeffs=[16.9952,0.0174535,-5.39619e-06,8.36867e-10,-5.22419e-14,-2750.07,-55.6309], Tmin=(1296.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.7317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][C]=CO(18753)',
    structure = SMILES('[CH2][C]=CO'),
    E0 = (186.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.23523,'amu*angstrom^2'), symmetry=1, barrier=(28.4004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2351,'amu*angstrom^2'), symmetry=1, barrier=(28.3973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24497,0.0260528,1.3484e-05,-5.00525e-08,2.54383e-11,22526.5,14.0801], Tmin=(100,'K'), Tmax=(898.827,'K')), NASAPolynomial(coeffs=[16.2027,-0.00210248,3.79693e-06,-8.3211e-10,5.63273e-14,18645.6,-59.4], Tmin=(898.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = 'C=[C]C(O)C=O(27744)',
    structure = SMILES('C=[C]C(O)C=O'),
    E0 = (-32.4957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,616.43],'cm^-1')),
        HinderedRotor(inertia=(0.478592,'amu*angstrom^2'), symmetry=1, barrier=(11.0038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0407669,'amu*angstrom^2'), symmetry=1, barrier=(10.9798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409213,'amu*angstrom^2'), symmetry=1, barrier=(10.9972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78798,0.043671,-3.37416e-05,1.31659e-08,-2.06728e-12,-3824.66,24.676], Tmin=(100,'K'), Tmax=(1501.82,'K')), NASAPolynomial(coeffs=[11.867,0.0168258,-6.92882e-06,1.26348e-09,-8.59241e-14,-6852.02,-28.0484], Tmin=(1501.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.4957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = 'OH(D)(132)',
    structure = SMILES('[OH]'),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92814e-05,-5.32177e-07,1.01951e-09,-3.85951e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.75,'K')), NASAPolynomial(coeffs=[3.07194,0.000604011,-1.39759e-08,-2.13452e-11,2.4807e-15,3579.39,4.57799], Tmin=(1145.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C=CC(C)[O](17702)',
    structure = SMILES('C=C=CC(C)[O]'),
    E0 = (181.611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,540,610,2055,3010,987.5,1337.5,450,1655,192.704,192.85],'cm^-1')),
        HinderedRotor(inertia=(0.654584,'amu*angstrom^2'), symmetry=1, barrier=(17.2625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.655299,'amu*angstrom^2'), symmetry=1, barrier=(17.2531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46721,0.0465386,-1.71976e-05,-1.05242e-08,7.08685e-12,21941.9,23.1054], Tmin=(100,'K'), Tmax=(1034.79,'K')), NASAPolynomial(coeffs=[12.7919,0.0216858,-8.6024e-06,1.61071e-09,-1.14458e-13,18585,-36.8125], Tmin=(1034.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C]C(O)[C](C)O(28170)',
    structure = SMILES('C=[C]C(O)[C](C)O'),
    E0 = (28.8905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622586,0.07768,-9.16449e-05,6.00393e-08,-1.60309e-11,3593.41,29.9795], Tmin=(100,'K'), Tmax=(906.89,'K')), NASAPolynomial(coeffs=[11.347,0.0303773,-1.34047e-05,2.52295e-09,-1.75248e-13,1648.27,-20.7112], Tmin=(906.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.8905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C=C(O)C(C)[O](13750)',
    structure = SMILES('[CH2]C=C(O)C(C)[O]'),
    E0 = (-58.0527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.452509,0.062859,-2.13948e-05,-2.82247e-08,1.83545e-11,-6840.59,28.0297], Tmin=(100,'K'), Tmax=(950.395,'K')), NASAPolynomial(coeffs=[20.6548,0.0148924,-4.1815e-06,7.30725e-10,-5.50816e-14,-12354.4,-77.2122], Tmin=(950.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.0527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=CC(O)C(C)[O](13756)',
    structure = SMILES('[CH]=CC(O)C(C)[O]'),
    E0 = (91.8778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,452.172,452.344],'cm^-1')),
        HinderedRotor(inertia=(0.0817265,'amu*angstrom^2'), symmetry=1, barrier=(11.8574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0816986,'amu*angstrom^2'), symmetry=1, barrier=(11.857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0816373,'amu*angstrom^2'), symmetry=1, barrier=(11.8576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0817622,'amu*angstrom^2'), symmetry=1, barrier=(11.8581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404337,0.0697194,-6.09578e-05,2.70887e-08,-4.76601e-12,11187.6,31.073], Tmin=(100,'K'), Tmax=(1373.76,'K')), NASAPolynomial(coeffs=[17.0363,0.0212916,-8.07954e-06,1.42755e-09,-9.61147e-14,6617.96,-54.4483], Tmin=(1373.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.8778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=C(O)C(C)O(28171)',
    structure = SMILES('[CH2][C]=C(O)C(C)O'),
    E0 = (-50.5717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.257934,0.0714421,-5.20774e-05,6.11214e-09,5.74369e-12,-5937.81,29.3935], Tmin=(100,'K'), Tmax=(957.924,'K')), NASAPolynomial(coeffs=[19.533,0.0161798,-5.04125e-06,8.66236e-10,-6.1307e-14,-10787.9,-68.8096], Tmin=(957.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.5717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)C(O)[C]=C(13923)',
    structure = SMILES('[CH2]C(O)C(O)[C]=C'),
    E0 = (63.8517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,400.172,400.173],'cm^-1')),
        HinderedRotor(inertia=(0.00105272,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0959038,'amu*angstrom^2'), symmetry=1, barrier=(10.8981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095903,'amu*angstrom^2'), symmetry=1, barrier=(10.8981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0959035,'amu*angstrom^2'), symmetry=1, barrier=(10.8981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0959048,'amu*angstrom^2'), symmetry=1, barrier=(10.8981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526406,0.0754046,-7.92648e-05,4.39261e-08,-9.75141e-12,7805.57,31.5259], Tmin=(100,'K'), Tmax=(1091.22,'K')), NASAPolynomial(coeffs=[14.1615,0.0254235,-1.05602e-05,1.95183e-09,-1.35022e-13,4829.81,-35.4457], Tmin=(1091.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.8517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=CC([O])C(C)[O](12747)',
    structure = SMILES('C=CC([O])C(C)[O]'),
    E0 = (75.1424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,402.302,402.95,403.877,403.888],'cm^-1')),
        HinderedRotor(inertia=(0.00103314,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131257,'amu*angstrom^2'), symmetry=1, barrier=(15.1716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132308,'amu*angstrom^2'), symmetry=1, barrier=(15.1637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4232,'J/mol'), sigma=(7.00504,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=661.03 K, Pc=27.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770968,0.0607713,-3.38901e-05,-4.14574e-10,4.59143e-12,9162.63,29.1558], Tmin=(100,'K'), Tmax=(1049.08,'K')), NASAPolynomial(coeffs=[15.3456,0.0245424,-9.74534e-06,1.81682e-09,-1.28467e-13,5040.25,-46.9299], Tmin=(1049.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.1424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC(O)[C](C)[O](13752)',
    structure = SMILES('C=CC(O)[C](C)[O]'),
    E0 = (21.4095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72099,0.0703653,-6.59376e-05,3.27333e-08,-6.59756e-12,2694.48,28.9508], Tmin=(100,'K'), Tmax=(1185.82,'K')), NASAPolynomial(coeffs=[13.4716,0.0273546,-1.15306e-05,2.14531e-09,-1.48778e-13,-329.451,-34.7362], Tmin=(1185.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.4095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C([O])C(O)C=C(13754)',
    structure = SMILES('[CH2]C([O])C(O)C=C'),
    E0 = (56.3707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,427.009,428.96,429.274],'cm^-1')),
        HinderedRotor(inertia=(0.000908045,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100005,'amu*angstrom^2'), symmetry=1, barrier=(13.1367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101948,'amu*angstrom^2'), symmetry=1, barrier=(13.1586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100751,'amu*angstrom^2'), symmetry=1, barrier=(13.1552,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.323032,0.071439,-6.43345e-05,2.93836e-08,-5.29416e-12,6920.11,31.5942], Tmin=(100,'K'), Tmax=(1346.61,'K')), NASAPolynomial(coeffs=[17.4302,0.0206232,-7.73003e-06,1.36026e-09,-9.1549e-14,2312.81,-56.0288], Tmin=(1346.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.3707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C]C([O])C(C)O(13755)',
    structure = SMILES('C=[C]C([O])C(C)O'),
    E0 = (82.6234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,386.026,386.617,388.566],'cm^-1')),
        HinderedRotor(inertia=(0.00111578,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110548,'amu*angstrom^2'), symmetry=1, barrier=(11.7223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111143,'amu*angstrom^2'), symmetry=1, barrier=(11.7451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110928,'amu*angstrom^2'), symmetry=1, barrier=(11.7536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.697401,0.0679222,-5.9482e-05,2.71022e-08,-4.97994e-12,10060.1,30.0851], Tmin=(100,'K'), Tmax=(1299.39,'K')), NASAPolynomial(coeffs=[14.5232,0.0253612,-1.03503e-05,1.89464e-09,-1.30057e-13,6467.03,-40.2374], Tmin=(1299.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.6234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=[C]C(O)C(C)O(28172)',
    structure = SMILES('[CH]=[C]C(O)C(C)O'),
    E0 = (99.3588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,304.655],'cm^-1')),
        HinderedRotor(inertia=(0.00181628,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168975,'amu*angstrom^2'), symmetry=1, barrier=(11.1291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168975,'amu*angstrom^2'), symmetry=1, barrier=(11.1291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168973,'amu*angstrom^2'), symmetry=1, barrier=(11.1291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168972,'amu*angstrom^2'), symmetry=1, barrier=(11.1291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.624956,0.0734979,-7.53082e-05,4.09814e-08,-8.98749e-12,12072.3,30.9418], Tmin=(100,'K'), Tmax=(1100.28,'K')), NASAPolynomial(coeffs=[13.5877,0.026372,-1.10611e-05,2.05317e-09,-1.4231e-13,9219.81,-32.8346], Tmin=(1100.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.3588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][CH]C(C)[O](17708)',
    structure = SMILES('C=[C][CH]C(C)[O]'),
    E0 = (372.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,399.072,399.599,399.695],'cm^-1')),
        HinderedRotor(inertia=(0.352754,'amu*angstrom^2'), symmetry=1, barrier=(39.9813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192794,'amu*angstrom^2'), symmetry=1, barrier=(21.8075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192533,'amu*angstrom^2'), symmetry=1, barrier=(21.8112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.28011,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08655,0.0583076,-4.75955e-05,1.99373e-08,-3.36764e-12,44892.6,22.6237], Tmin=(100,'K'), Tmax=(1404.68,'K')), NASAPolynomial(coeffs=[13.7617,0.0222135,-9.05218e-06,1.64451e-09,-1.11963e-13,41331.6,-42.834], Tmin=(1404.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C=[C]C([O])C(C)[O](13760)',
    structure = SMILES('C=[C]C([O])C(C)[O]'),
    E0 = (312.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2950,3100,1380,975,1025,1650,409.73,410.078,410.529,410.772],'cm^-1')),
        HinderedRotor(inertia=(0.00100242,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114223,'amu*angstrom^2'), symmetry=1, barrier=(13.5956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113882,'amu*angstrom^2'), symmetry=1, barrier=(13.5949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803562,0.0645657,-5.5657e-05,2.45519e-08,-4.34299e-12,37763.1,29.5141], Tmin=(100,'K'), Tmax=(1350.85,'K')), NASAPolynomial(coeffs=[14.8711,0.0229105,-9.40281e-06,1.72473e-09,-1.18418e-13,33962.4,-42.5842], Tmin=(1350.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)[CH][O](28102)',
    structure = SMILES('C=[C]C(O)[CH][O]'),
    E0 = (294.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,276.259,279.06,2001.85],'cm^-1')),
        HinderedRotor(inertia=(0.00227246,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191772,'amu*angstrom^2'), symmetry=1, barrier=(10.9158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189866,'amu*angstrom^2'), symmetry=1, barrier=(10.8672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40118,0.0638678,-0.000101269,9.07849e-08,-3.19727e-11,35530,26.2595], Tmin=(100,'K'), Tmax=(835.883,'K')), NASAPolynomial(coeffs=[6.37686,0.0274382,-1.32505e-05,2.52384e-09,-1.73332e-13,35139,5.78378], Tmin=(835.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S) + radical(CCsJOH)"""),
)

species(
    label = 'C=[C]C(O)[C](C)[O](28173)',
    structure = SMILES('C=[C]C(O)[C](C)[O]'),
    E0 = (259.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,206.25,206.282,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00396257,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569591,'amu*angstrom^2'), symmetry=1, barrier=(17.1921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569342,'amu*angstrom^2'), symmetry=1, barrier=(17.1923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569547,'amu*angstrom^2'), symmetry=1, barrier=(17.1924,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.804348,0.0734637,-8.49639e-05,5.39971e-08,-1.39909e-11,31293.1,29.1353], Tmin=(100,'K'), Tmax=(932.629,'K')), NASAPolynomial(coeffs=[11.2604,0.0286184,-1.28373e-05,2.4396e-09,-1.70592e-13,29342.8,-20.5802], Tmin=(932.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2][C]=C(O)C(C)[O](28174)',
    structure = SMILES('[CH2][C]=C(O)C(C)[O]'),
    E0 = (179.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,422,422.15],'cm^-1')),
        HinderedRotor(inertia=(0.128499,'amu*angstrom^2'), symmetry=1, barrier=(16.2534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128506,'amu*angstrom^2'), symmetry=1, barrier=(16.2526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128573,'amu*angstrom^2'), symmetry=1, barrier=(16.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128435,'amu*angstrom^2'), symmetry=1, barrier=(16.2514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.418081,0.0674797,-4.62809e-05,1.21988e-09,7.28759e-12,21762.8,28.6268], Tmin=(100,'K'), Tmax=(959.405,'K')), NASAPolynomial(coeffs=[19.4749,0.0143755,-4.44916e-06,7.77322e-10,-5.62065e-14,16893.6,-68.8411], Tmin=(959.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([O])C(O)[C]=C(14793)',
    structure = SMILES('[CH2]C([O])C(O)[C]=C'),
    E0 = (294.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,427.651,427.885,428.207],'cm^-1')),
        HinderedRotor(inertia=(0.000922752,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916945,'amu*angstrom^2'), symmetry=1, barrier=(11.9114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916249,'amu*angstrom^2'), symmetry=1, barrier=(11.9094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0917068,'amu*angstrom^2'), symmetry=1, barrier=(11.9019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671625,0.071628,-7.41515e-05,3.99426e-08,-8.59586e-12,35506.8,30.812], Tmin=(100,'K'), Tmax=(1125.09,'K')), NASAPolynomial(coeffs=[14.2408,0.0233856,-9.83334e-06,1.83106e-09,-1.27278e-13,32453.5,-36.2509], Tmin=(1125.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(O)C(C)[O](28175)',
    structure = SMILES('[CH]=[C]C(O)C(C)[O]'),
    E0 = (329.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,336.784,340.604],'cm^-1')),
        HinderedRotor(inertia=(0.00146322,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14965,'amu*angstrom^2'), symmetry=1, barrier=(12.1507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149663,'amu*angstrom^2'), symmetry=1, barrier=(12.1101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151104,'amu*angstrom^2'), symmetry=1, barrier=(12.1176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769319,0.0697327,-7.02388e-05,3.70571e-08,-7.85705e-12,39773.6,30.2308], Tmin=(100,'K'), Tmax=(1137.17,'K')), NASAPolynomial(coeffs=[13.6837,0.0243061,-1.03182e-05,1.92863e-09,-1.34254e-13,36836.4,-33.7339], Tmin=(1137.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC(O)C(C)=O(13763)',
    structure = SMILES('C=CC(O)C(C)=O'),
    E0 = (-325.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0618,0.0556095,-2.99883e-05,1.34976e-09,2.75317e-12,-39000.2,28.3987], Tmin=(100,'K'), Tmax=(1110.77,'K')), NASAPolynomial(coeffs=[13.3828,0.0262972,-1.07381e-05,2.00002e-09,-1.39935e-13,-42666.2,-36.5188], Tmin=(1110.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-325.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C(O)C(C)O(28176)',
    structure = SMILES('C=C=C(O)C(C)O'),
    E0 = (-263.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.33177,0.0689642,-4.48301e-05,-1.44485e-09,8.41019e-12,-31526,27.4279], Tmin=(100,'K'), Tmax=(956.809,'K')), NASAPolynomial(coeffs=[19.6269,0.0158855,-4.8648e-06,8.41474e-10,-6.03743e-14,-36481.1,-71.4067], Tmin=(956.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C=[C]C(O)C[O](13587)',
    structure = SMILES('C=[C]C(O)C[O]'),
    E0 = (114.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,217.895,221.788,4000],'cm^-1')),
        HinderedRotor(inertia=(0.335038,'amu*angstrom^2'), symmetry=1, barrier=(11.2479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329155,'amu*angstrom^2'), symmetry=1, barrier=(11.2154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141791,'amu*angstrom^2'), symmetry=1, barrier=(30.0625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4073.83,'J/mol'), sigma=(6.66081,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=636.32 K, Pc=31.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25621,0.0458497,-2.15226e-05,-4.3278e-08,5.33123e-11,13813.2,23.7683], Tmin=(100,'K'), Tmax=(475.481,'K')), NASAPolynomial(coeffs=[5.50376,0.031112,-1.47234e-05,2.84319e-09,-1.99498e-13,13362.1,9.01943], Tmin=(475.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S)"""),
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
    label = 'C=C1OC(C)C1O(28158)',
    structure = SMILES('C=C1OC(C)C1O'),
    E0 = (-264.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05332,0.0415413,4.71464e-05,-1.07902e-07,4.99267e-11,-31680.2,19.967], Tmin=(100,'K'), Tmax=(909.457,'K')), NASAPolynomial(coeffs=[22.7001,0.00766501,1.86427e-06,-5.57935e-10,3.564e-14,-38154,-96.356], Tmin=(909.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'C=[C][CH]C(C)OO(28177)',
    structure = SMILES('C=[C][CH]C(C)OO'),
    E0 = (213.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,296.18,296.184],'cm^-1')),
        HinderedRotor(inertia=(0.350241,'amu*angstrom^2'), symmetry=1, barrier=(21.8022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.350242,'amu*angstrom^2'), symmetry=1, barrier=(21.8022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.350244,'amu*angstrom^2'), symmetry=1, barrier=(21.8023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01579,'amu*angstrom^2'), symmetry=1, barrier=(63.2355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.604412,'amu*angstrom^2'), symmetry=1, barrier=(37.6242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615648,0.0749278,-6.97765e-05,3.42876e-08,-6.91992e-12,25798.5,25.2663], Tmin=(100,'K'), Tmax=(1172.42,'K')), NASAPolynomial(coeffs=[13.3261,0.031563,-1.42953e-05,2.73962e-09,-1.92819e-13,22818.1,-38.0761], Tmin=(1172.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJCO)"""),
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
    label = 'C=[C]C(O)[CH]C(27811)',
    structure = SMILES('C=[C]C(O)[CH]C'),
    E0 = (224.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,398.902,401.27],'cm^-1')),
        HinderedRotor(inertia=(0.00105341,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00105824,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867638,'amu*angstrom^2'), symmetry=1, barrier=(9.88923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0864274,'amu*angstrom^2'), symmetry=1, barrier=(9.87125,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34412,0.0504701,-2.70394e-05,2.37526e-10,3.18803e-12,27119.8,28.0942], Tmin=(100,'K'), Tmax=(1064.42,'K')), NASAPolynomial(coeffs=[12.3302,0.0234925,-9.18409e-06,1.68217e-09,-1.17125e-13,23970.5,-29.4006], Tmin=(1064.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[C]=C(584)',
    structure = SMILES('[C]=C'),
    E0 = (600.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94093,-0.00117598,1.80376e-05,-2.01208e-08,6.96659e-12,72197.9,5.25681], Tmin=(100,'K'), Tmax=(976.125,'K')), NASAPolynomial(coeffs=[3.93016,0.00536132,-1.98619e-06,3.69549e-10,-2.66221e-14,71890.7,3.724], Tmin=(976.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'CC([O])[CH]O(3785)',
    structure = SMILES('CC([O])[CH]O'),
    E0 = (-42.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,322.464,322.474],'cm^-1')),
        HinderedRotor(inertia=(0.00162117,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132686,'amu*angstrom^2'), symmetry=1, barrier=(9.78974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132702,'amu*angstrom^2'), symmetry=1, barrier=(9.78979,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68986,0.0521771,-5.85792e-05,3.66864e-08,-9.35822e-12,-4989.01,20.0755], Tmin=(100,'K'), Tmax=(948.754,'K')), NASAPolynomial(coeffs=[9.18448,0.0205795,-8.62305e-06,1.58361e-09,-1.08566e-13,-6411.13,-15.6876], Tmin=(948.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCsJOH)"""),
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
    E0 = (82.6234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (195.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (191.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (235.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (82.6234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (134.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (145.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (210.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (248.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (265.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (197.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (157.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (221.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (224.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (224.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (126.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (188.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (132.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (400.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (531.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (336.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (430.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (471.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (391.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (506.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (541.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (171.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (160.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (533.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (177.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (90.9078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (320.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (631.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (592.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['CC=O(606)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', 'C=[C]C(O)C(C)=O(27722)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C=C=C(O)C(C)[O](28168)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC(O)C(C)[O](28169)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CC=O(606)', '[CH2][C]=CO(18753)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(73.8577,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 70.1 to 73.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH3](11)', 'C=[C]C(O)C=O(27744)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-CsH_O;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C[CH][O](605)', 'C=C=CO(12571)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00401797,'m^3/(mol*s)'), n=2.41733, Ea=(22.1495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(D)(132)', 'C=C=CC(C)[O](17702)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(986500,'cm^3/(mol*s)'), n=2.037, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;OJ_pri] for rate rule [Cds-CsH_Ca;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['C=[C]C(O)[C](C)O(28170)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['[CH2]C=C(O)C(C)[O](13750)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(O)C(C)[O](13756)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['[CH2][C]=C(O)C(C)O(28171)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['[CH2]C(O)C(O)[C]=C(13923)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['C=CC([O])C(C)[O](12747)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['C=CC(O)[C](C)[O](13752)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['[CH2]C([O])C(O)C=C(13754)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C([O])C(C)O(13755)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C(O)C(C)O(28172)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OH(D)(132)', 'C=[C][CH]C(C)[O](17708)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', 'C=[C]C([O])C(C)[O](13760)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.00518e+06,'m^3/(mol*s)'), n=0.282325, Ea=(7.09479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH][O](605)', '[CH2][C]=CO(18753)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH3](11)', 'C=[C]C(O)[CH][O](28102)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', 'C=[C]C(O)[C](C)[O](28173)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH2][C]=C(O)C(C)[O](28174)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH]=[C]C(O)C(C)[O](28175)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['C=CC(O)C(C)=O(13763)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['C=C=C(O)C(C)O(28176)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(14)', 'C=[C]C(O)C[O](13587)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['[CH2]C(=CO)C(C)[O](13701)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]C(O)C(C)[O](13753)'],
    products = ['C=C1OC(C)C1O(28158)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C][CH]C(C)OO(28177)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(T)(63)', 'C=[C]C(O)[CH]C(27811)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[C]=C(584)', 'CC([O])[CH]O(3785)'],
    products = ['C=[C]C(O)C(C)[O](13753)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

network(
    label = '5099',
    isomers = [
        'C=[C]C(O)C(C)[O](13753)',
    ],
    reactants = [
        ('CC=O(606)', 'C=C=CO(12571)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5099',
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

