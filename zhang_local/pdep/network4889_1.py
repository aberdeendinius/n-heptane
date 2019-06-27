species(
    label = '[O]C=C=CCC=C[O](22588)',
    structure = SMILES('[O]C=C=CCC=C[O]'),
    E0 = (94.6684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24294,'amu*angstrom^2'), symmetry=1, barrier=(28.5775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24272,'amu*angstrom^2'), symmetry=1, barrier=(28.5727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371633,0.0609434,-8.25709e-06,-4.93002e-08,2.78157e-11,11534,29.2454], Tmin=(100,'K'), Tmax=(935.092,'K')), NASAPolynomial(coeffs=[23.8145,0.00848651,-8.25121e-07,9.47677e-11,-1.27022e-14,5058.87,-93.4597], Tmin=(935.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.6684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = 'C#CC=O(21959)',
    structure = SMILES('C#CC=O'),
    E0 = (84.2941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,346.083,346.148,1254.41,1569.49,1569.49],'cm^-1')),
        HinderedRotor(inertia=(0.289847,'amu*angstrom^2'), symmetry=1, barrier=(24.6472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3153.34,'J/mol'), sigma=(5.09602,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.54 K, Pc=54.07 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00307,0.020875,-1.61807e-05,6.30357e-09,-1.00066e-12,10174.9,9.76778], Tmin=(100,'K'), Tmax=(1462.04,'K')), NASAPolynomial(coeffs=[7.33456,0.00902438,-4.02236e-06,7.59526e-10,-5.26541e-14,8908.32,-12.7744], Tmin=(1462.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.2941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = '[O]C=[C]C1CC=CO1(25557)',
    structure = SMILES('[O]C=[C]C1CC=CO1'),
    E0 = (109.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09534,0.0336498,8.0287e-05,-1.52632e-07,6.86225e-11,13270.3,23.3441], Tmin=(100,'K'), Tmax=(905.937,'K')), NASAPolynomial(coeffs=[27.2296,-0.00194685,7.10682e-06,-1.55563e-09,1.01919e-13,5260.68,-118.229], Tmin=(905.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'O=CC1C[CH][C]=CO1(25456)',
    structure = SMILES('O=CC1C[CH][C]=CO1'),
    E0 = (107.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73278,0.0234154,8.91851e-05,-1.43262e-07,5.98258e-11,13087.6,21.3915], Tmin=(100,'K'), Tmax=(927.908,'K')), NASAPolynomial(coeffs=[20.3471,0.0110385,-5.13573e-07,3.24887e-12,-9.06152e-15,6711.46,-82.762], Tmin=(927.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(3,4-Dihydro-2H-pyran) + radical(Allyl_S) + radical(Cds_S)"""),
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
    label = '[O]C=C=CCC=C=O(25685)',
    structure = SMILES('[O]C=C=CCC=C=O'),
    E0 = (75.5147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0445,'amu*angstrom^2'), symmetry=1, barrier=(24.0151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04177,'amu*angstrom^2'), symmetry=1, barrier=(23.9524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424297,0.0732837,-7.67944e-05,4.04966e-08,-8.36189e-12,9215.63,27.0498], Tmin=(100,'K'), Tmax=(1185.12,'K')), NASAPolynomial(coeffs=[16.7481,0.0181877,-7.05985e-06,1.26875e-09,-8.68146e-14,5346.49,-54.4757], Tmin=(1185.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.5147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CCC=C=C=O(25686)',
    structure = SMILES('[O]C=CCC=C=C=O'),
    E0 = (138.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0926077,'amu*angstrom^2'), symmetry=1, barrier=(2.12923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0924983,'amu*angstrom^2'), symmetry=1, barrier=(2.12672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65228,0.0373643,1.01286e-05,-4.73659e-08,2.27535e-11,16787.1,11.8507], Tmin=(100,'K'), Tmax=(952.59,'K')), NASAPolynomial(coeffs=[16.6143,0.00994566,-2.45196e-06,4.58895e-10,-3.83125e-14,12330,-68.0382], Tmin=(952.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CCC=[C]O(25687)',
    structure = SMILES('[O]C=C=CCC=[C]O'),
    E0 = (192.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0116,'amu*angstrom^2'), symmetry=1, barrier=(23.2586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00756,'amu*angstrom^2'), symmetry=1, barrier=(23.1659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00711,'amu*angstrom^2'), symmetry=1, barrier=(23.1554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169504,0.0696731,-3.85264e-05,-1.67955e-08,1.63885e-11,23357.9,31.6896], Tmin=(100,'K'), Tmax=(926.396,'K')), NASAPolynomial(coeffs=[22.9819,0.008834,-9.96421e-07,7.94876e-11,-7.81948e-15,17515.2,-85.3453], Tmin=(926.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CC[CH]C#CO(25688)',
    structure = SMILES('[O]C=CC[CH]C#CO'),
    E0 = (167.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2100,2250,500,550,2995,3025,975,1000,1300,1375,400,500,1630,1680,405.425,405.955,406.124,406.135],'cm^-1')),
        HinderedRotor(inertia=(0.778435,'amu*angstrom^2'), symmetry=1, barrier=(91.2472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165118,'amu*angstrom^2'), symmetry=1, barrier=(19.2726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165223,'amu*angstrom^2'), symmetry=1, barrier=(19.276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781639,'amu*angstrom^2'), symmetry=1, barrier=(91.2481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.175706,0.0738997,-7.02906e-05,3.34106e-08,-6.02602e-12,20277.4,32.3277], Tmin=(100,'K'), Tmax=(1542.19,'K')), NASAPolynomial(coeffs=[19.1959,0.0141699,-2.9688e-06,3.20059e-10,-1.52788e-14,15430.4,-65.8637], Tmin=(1542.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-CtH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=COJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[O]C=C=CC[C]=CO(25689)',
    structure = SMILES('[O]C=C=CC[C]=CO'),
    E0 = (191.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11134,'amu*angstrom^2'), symmetry=1, barrier=(25.5519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11128,'amu*angstrom^2'), symmetry=1, barrier=(25.5506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11143,'amu*angstrom^2'), symmetry=1, barrier=(25.554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0688376,0.0724563,-3.61978e-05,-2.66638e-08,2.18221e-11,23140.1,29.8906], Tmin=(100,'K'), Tmax=(916.915,'K')), NASAPolynomial(coeffs=[25.9182,0.00423282,1.55867e-06,-4.19886e-10,2.63005e-14,16476.8,-103.576], Tmin=(916.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C=C[CH]C=CO(25690)',
    structure = SMILES('[O]C=C=C[CH]C=CO'),
    E0 = (54.9308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277228,0.0555203,2.61182e-05,-1.00683e-07,5.07216e-11,6765.37,27.8461], Tmin=(100,'K'), Tmax=(908.931,'K')), NASAPolynomial(coeffs=[29.3612,-0.0024519,6.23615e-06,-1.34677e-09,8.77881e-14,-1414.07,-125.601], Tmin=(908.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.9308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CCC#C[CH]O(25691)',
    structure = SMILES('[O]C=CCC#C[CH]O'),
    E0 = (175.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2100,2250,500,550,2995,3025,975,1000,1300,1375,400,500,1630,1680,410.905,410.914,410.952,411.014],'cm^-1')),
        HinderedRotor(inertia=(0.152275,'amu*angstrom^2'), symmetry=1, barrier=(18.2328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122976,'amu*angstrom^2'), symmetry=1, barrier=(14.7324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122951,'amu*angstrom^2'), symmetry=1, barrier=(14.7324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622372,'amu*angstrom^2'), symmetry=1, barrier=(74.5674,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0487799,0.0755545,-7.75633e-05,4.04255e-08,-8.09809e-12,21293.5,32.1655], Tmin=(100,'K'), Tmax=(1328.8,'K')), NASAPolynomial(coeffs=[18.2103,0.015618,-3.96016e-06,5.15919e-10,-2.83917e-14,16931.8,-58.8667], Tmin=(1328.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=COJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C=C=[C]CC=CO(25692)',
    structure = SMILES('[O]C=C=[C]CC=CO'),
    E0 = (191.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11134,'amu*angstrom^2'), symmetry=1, barrier=(25.5519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11128,'amu*angstrom^2'), symmetry=1, barrier=(25.5506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11143,'amu*angstrom^2'), symmetry=1, barrier=(25.554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0688376,0.0724563,-3.61978e-05,-2.66638e-08,2.18221e-11,23140.1,29.8906], Tmin=(100,'K'), Tmax=(916.915,'K')), NASAPolynomial(coeffs=[25.9182,0.00423282,1.55867e-06,-4.19886e-10,2.63005e-14,16476.8,-103.576], Tmin=(916.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C[CH]C=C=CO(25693)',
    structure = SMILES('[O]C=C[CH]C=C=CO'),
    E0 = (54.9308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277228,0.0555203,2.61182e-05,-1.00683e-07,5.07216e-11,6765.37,27.8461], Tmin=(100,'K'), Tmax=(908.931,'K')), NASAPolynomial(coeffs=[29.3612,-0.0024519,6.23615e-06,-1.34677e-09,8.77881e-14,-1414.07,-125.601], Tmin=(908.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.9308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = '[O]C=[C]CC=C=CO(25694)',
    structure = SMILES('[O]C=[C]CC=C=CO'),
    E0 = (191.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11134,'amu*angstrom^2'), symmetry=1, barrier=(25.5519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11128,'amu*angstrom^2'), symmetry=1, barrier=(25.5506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11143,'amu*angstrom^2'), symmetry=1, barrier=(25.554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0688376,0.0724563,-3.61978e-05,-2.66638e-08,2.18221e-11,23140.1,29.8906], Tmin=(100,'K'), Tmax=(916.915,'K')), NASAPolynomial(coeffs=[25.9182,0.00423282,1.55867e-06,-4.19886e-10,2.63005e-14,16476.8,-103.576], Tmin=(916.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'O=C=[C][CH]CC=CO(25695)',
    structure = SMILES('O=C=[C][CH]CC=CO'),
    E0 = (119.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05798,'amu*angstrom^2'), symmetry=1, barrier=(24.3251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05848,'amu*angstrom^2'), symmetry=1, barrier=(24.3365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.058,'amu*angstrom^2'), symmetry=1, barrier=(24.3256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0596,'amu*angstrom^2'), symmetry=1, barrier=(24.3624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.377521,0.0817048,-8.38588e-05,4.18731e-08,-7.97521e-12,14521.9,31.8005], Tmin=(100,'K'), Tmax=(1386.82,'K')), NASAPolynomial(coeffs=[21.854,0.0124194,-3.33443e-06,4.7915e-10,-2.92203e-14,8852.24,-80.9334], Tmin=(1386.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCCJ=C=O) + radical(CCJC(C)=C=O)"""),
)

species(
    label = 'O=[C][CH]CC=C=CO(25696)',
    structure = SMILES('O=[C][CH]CC=C=CO'),
    E0 = (136.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,540,610,2055,1855,455,950,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.871068,'amu*angstrom^2'), symmetry=1, barrier=(20.0276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871055,'amu*angstrom^2'), symmetry=1, barrier=(20.0273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.870665,'amu*angstrom^2'), symmetry=1, barrier=(20.0183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.870759,'amu*angstrom^2'), symmetry=1, barrier=(20.0205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.511796,0.0825141,-8.75057e-05,4.48906e-08,-8.64994e-12,16554.8,31.8493], Tmin=(100,'K'), Tmax=(1435.57,'K')), NASAPolynomial(coeffs=[21.9412,0.00999171,-1.3208e-06,3.38457e-11,3.41689e-15,11134.6,-81.0167], Tmin=(1435.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C[O](602)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (221.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27415,0.00479611,3.742e-05,-6.11894e-08,2.65903e-11,26726.9,9.63858], Tmin=(100,'K'), Tmax=(905.806,'K')), NASAPolynomial(coeffs=[11.9892,-0.00434473,3.96329e-06,-8.00891e-10,5.23184e-14,23944.2,-38.1893], Tmin=(905.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = 'C=C[C]=C[O](22438)',
    structure = SMILES('C=C[C]=C[O]'),
    E0 = (225.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61747,'amu*angstrom^2'), symmetry=1, barrier=(37.1889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98077,0.0322088,7.53766e-06,-4.54875e-08,2.38742e-11,27216,16.1107], Tmin=(100,'K'), Tmax=(899.949,'K')), NASAPolynomial(coeffs=[16.1069,0.00249504,1.93925e-06,-5.05324e-10,3.47549e-14,23334.1,-57.9915], Tmin=(899.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ)"""),
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
    label = '[CH]=C=C[O](8556)',
    structure = SMILES('[CH]=C=C[O]'),
    E0 = (269.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7401,0.0176211,1.62543e-05,-4.58526e-08,2.291e-11,32513.4,12.64], Tmin=(100,'K'), Tmax=(883.628,'K')), NASAPolynomial(coeffs=[13.5244,-0.00323064,4.17673e-06,-9.22668e-10,6.45049e-14,29515.7,-44.2318], Tmin=(883.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=C[CH]C=C[O](25697)',
    structure = SMILES('[O]C=C=C[CH]C=C[O]'),
    E0 = (196.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62895,'amu*angstrom^2'), symmetry=1, barrier=(37.4527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62046,'amu*angstrom^2'), symmetry=1, barrier=(37.2576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695524,0.0484862,2.96475e-05,-9.44386e-08,4.58555e-11,23762.2,27.7539], Tmin=(100,'K'), Tmax=(918.912,'K')), NASAPolynomial(coeffs=[26.001,0.00141529,3.50953e-06,-7.67763e-10,4.61873e-14,16448.1,-106.682], Tmin=(918.912,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[O]C=[C]CC=C=C[O](25698)',
    structure = SMILES('[O]C=[C]CC=C=C[O]'),
    E0 = (332.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2191,'amu*angstrom^2'), symmetry=1, barrier=(28.0296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22259,'amu*angstrom^2'), symmetry=1, barrier=(28.1097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342925,0.0654961,-3.29055e-05,-2.01592e-08,1.68742e-11,40137.1,29.822], Tmin=(100,'K'), Tmax=(936.243,'K')), NASAPolynomial(coeffs=[22.606,0.00801836,-1.12086e-06,1.48001e-10,-1.43778e-14,34318.8,-84.9277], Tmin=(936.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=[C]CC=C[O](25699)',
    structure = SMILES('[O]C=C=[C]CC=C[O]'),
    E0 = (332.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2191,'amu*angstrom^2'), symmetry=1, barrier=(28.0296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22259,'amu*angstrom^2'), symmetry=1, barrier=(28.1097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342925,0.0654961,-3.29055e-05,-2.01592e-08,1.68742e-11,40137.1,29.822], Tmin=(100,'K'), Tmax=(936.243,'K')), NASAPolynomial(coeffs=[22.606,0.00801836,-1.12086e-06,1.48001e-10,-1.43778e-14,34318.8,-84.9277], Tmin=(936.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CC[CH][C]=O(25700)',
    structure = SMILES('[O]C=C=CC[CH][C]=O'),
    E0 = (277.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,298.399,298.402,298.403],'cm^-1')),
        HinderedRotor(inertia=(0.32459,'amu*angstrom^2'), symmetry=1, barrier=(20.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.324574,'amu*angstrom^2'), symmetry=1, barrier=(20.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.324584,'amu*angstrom^2'), symmetry=1, barrier=(20.5094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.157001,0.0724768,-7.31833e-05,3.66007e-08,-7.03051e-12,33540.7,30.861], Tmin=(100,'K'), Tmax=(1360.62,'K')), NASAPolynomial(coeffs=[19.0109,0.0132833,-3.77394e-06,5.57691e-10,-3.41675e-14,28758.7,-64.6227], Tmin=(1360.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCJCHO) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=CC[CH][C]=C=O(25701)',
    structure = SMILES('[O]C=CC[CH][C]=C=O'),
    E0 = (260.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,322.398,322.429,322.454,322.481],'cm^-1')),
        HinderedRotor(inertia=(0.37298,'amu*angstrom^2'), symmetry=1, barrier=(27.5176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37291,'amu*angstrom^2'), symmetry=1, barrier=(27.5149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373178,'amu*angstrom^2'), symmetry=1, barrier=(27.5166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315252,0.0714184,-6.8814e-05,3.28397e-08,-6.11181e-12,31506.7,30.724], Tmin=(100,'K'), Tmax=(1314.06,'K')), NASAPolynomial(coeffs=[18.3145,0.0166288,-6.27171e-06,1.10997e-09,-7.52346e-14,26776.3,-61.0279], Tmin=(1314.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCJC(C)=C=O) + radical(CCCJ=C=O) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CCC1[CH]O1(25702)',
    structure = SMILES('[O]C=C=CCC1[CH]O1'),
    E0 = (233.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26742,0.0837699,-8.4178e-05,4.11724e-08,-7.33892e-12,28332.5,32.979], Tmin=(100,'K'), Tmax=(1681.56,'K')), NASAPolynomial(coeffs=[20.1475,0.0104271,1.22928e-06,-6.10503e-10,5.08704e-14,24297.7,-72.0479], Tmin=(1681.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Ethylene_oxide) + radical(C=COJ) + radical(CCsJO)"""),
)

species(
    label = '[O]C=CCC=C1[CH]O1(25703)',
    structure = SMILES('[O]C=CCC=C1[CH]O1'),
    E0 = (120.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.758377,0.0401444,6.76675e-05,-1.40207e-07,6.31756e-11,14680.6,26.439], Tmin=(100,'K'), Tmax=(924.338,'K')), NASAPolynomial(coeffs=[29.0899,-0.00196256,5.37116e-06,-1.06372e-09,6.12168e-14,6004.22,-126.616], Tmin=(924.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(methyleneoxirane) + radical(C=COJ) + radical(C=CCJO)"""),
)

species(
    label = '[O]C=CCC1[C]=CO1(25704)',
    structure = SMILES('[O]C=CCC1[C]=CO1'),
    E0 = (213.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5358,0.0456415,5.52791e-05,-1.29128e-07,5.95499e-11,25799.5,26.9175], Tmin=(100,'K'), Tmax=(925.362,'K')), NASAPolynomial(coeffs=[29.9136,-0.00246227,5.38206e-06,-1.05575e-09,6.06286e-14,16985,-130.784], Tmin=(925.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C1[CH]CC=CO1(25504)',
    structure = SMILES('[O]C=C1[CH]CC=CO1'),
    E0 = (76.8356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65778,0.0201669,0.000110247,-1.65208e-07,6.59866e-11,9354.01,23.8878], Tmin=(100,'K'), Tmax=(956.226,'K')), NASAPolynomial(coeffs=[22.2043,0.013466,-3.55458e-06,8.01798e-10,-7.50291e-14,1801.51,-93.2624], Tmin=(956.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.8356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclohexane) + radical(CCJCO) + radical(C=COJ)"""),
)

species(
    label = '[O]C1[C]=CCC=CO1(25636)',
    structure = SMILES('[O]C1[C]=CCC=CO1'),
    E0 = (216.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41007,-0.000561887,0.000223504,-3.36522e-07,1.43141e-10,26139.4,25.8232], Tmin=(100,'K'), Tmax=(903.862,'K')), NASAPolynomial(coeffs=[42.4317,-0.0312492,2.40859e-05,-4.78676e-09,3.14064e-13,12561.8,-202.023], Tmin=(903.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[O]C1[CH]CC=C=CO1(25705)',
    structure = SMILES('[O]C1[CH]CC=C=CO1'),
    E0 = (276.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24667,0.0433058,2.1149e-05,-6.07437e-08,2.64695e-11,33380.4,21.796], Tmin=(100,'K'), Tmax=(986.367,'K')), NASAPolynomial(coeffs=[16.9376,0.0210989,-8.0755e-06,1.58588e-09,-1.19714e-13,28269.8,-63.9033], Tmin=(986.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'O=C=CCC=C=CO(25706)',
    structure = SMILES('O=C=CCC=C=CO'),
    E0 = (-65.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0082431,0.0805379,-8.14056e-05,3.62432e-08,-4.66511e-12,-7780.52,27.1904], Tmin=(100,'K'), Tmax=(960.641,'K')), NASAPolynomial(coeffs=[19.6382,0.0150954,-4.77046e-06,7.91409e-10,-5.35481e-14,-12310.2,-70.7335], Tmin=(960.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O=C=C=CCC=CO(25707)',
    structure = SMILES('O=C=C=CCC=CO'),
    E0 = (-2.69765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23642,0.0443849,6.55999e-06,-5.33988e-08,2.74403e-11,-209.807,11.9333], Tmin=(100,'K'), Tmax=(929.65,'K')), NASAPolynomial(coeffs=[19.8995,0.00620394,2.03163e-07,-1.03377e-10,1.90958e-15,-5499.98,-86.5333], Tmin=(929.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.69765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[O][CH]CC=C[C]=C[O](25708)',
    structure = SMILES('[O][CH]CC=C[C]=C[O]'),
    E0 = (413.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,230.039,230.236,230.241,230.323,230.426],'cm^-1')),
        HinderedRotor(inertia=(0.506805,'amu*angstrom^2'), symmetry=1, barrier=(19.077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506917,'amu*angstrom^2'), symmetry=1, barrier=(19.0806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02441,'amu*angstrom^2'), symmetry=1, barrier=(76.1785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.1637,0.0840579,-9.95654e-05,6.18469e-08,-1.51912e-11,49821.9,30.4665], Tmin=(100,'K'), Tmax=(997.234,'K')), NASAPolynomial(coeffs=[15.0722,0.0242573,-9.61371e-06,1.71152e-09,-1.1534e-13,46848.5,-41.4168], Tmin=(997.234,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C=[C]C[CH]C=C[O](25709)',
    structure = SMILES('[O]C=[C]C[CH]C=C[O]'),
    E0 = (307.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,343.733,343.734,343.734,343.734,343.734],'cm^-1')),
        HinderedRotor(inertia=(0.301502,'amu*angstrom^2'), symmetry=1, barrier=(25.2791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301502,'amu*angstrom^2'), symmetry=1, barrier=(25.2791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301503,'amu*angstrom^2'), symmetry=1, barrier=(25.2791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439369,0.0607057,-1.14695e-05,-4.28663e-08,2.47041e-11,37187.5,29.7268], Tmin=(100,'K'), Tmax=(940.641,'K')), NASAPolynomial(coeffs=[22.5689,0.0107131,-2.09081e-06,3.40932e-10,-2.93994e-14,31072.8,-86.0543], Tmin=(940.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(C=COJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C[C]CC=C[O](25710)',
    structure = SMILES('[O]C=C[C]CC=C[O]'),
    E0 = (377.255,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18828,'amu*angstrom^2'), symmetry=1, barrier=(27.321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18837,'amu*angstrom^2'), symmetry=1, barrier=(27.3229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18807,'amu*angstrom^2'), symmetry=1, barrier=(27.3162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.221388,0.064385,-1.51117e-05,-4.40913e-08,2.64103e-11,45526.6,31.0479], Tmin=(100,'K'), Tmax=(933.793,'K')), NASAPolynomial(coeffs=[24.5564,0.00771451,-4.94938e-07,2.92032e-11,-7.91346e-15,38907.8,-95.7929], Tmin=(933.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCJ2_triplet) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=C[CH]CC=C[O](25711)',
    structure = SMILES('[O][C]=C[CH]CC=C[O]'),
    E0 = (309.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,474.863,474.863,474.879,474.979,475.047],'cm^-1')),
        HinderedRotor(inertia=(0.135507,'amu*angstrom^2'), symmetry=1, barrier=(21.6944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135593,'amu*angstrom^2'), symmetry=1, barrier=(21.694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135512,'amu*angstrom^2'), symmetry=1, barrier=(21.6943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675182,0.0579433,-1.38191e-05,-3.30661e-08,1.93519e-11,37405.5,31.5354], Tmin=(100,'K'), Tmax=(955.76,'K')), NASAPolynomial(coeffs=[19.6833,0.0152296,-4.59774e-06,8.29037e-10,-6.25917e-14,32089.5,-68.1096], Tmin=(955.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(C=COJ) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[O]C=[C]CC[C]=C[O](25712)',
    structure = SMILES('[O]C=[C]CC[C]=C[O]'),
    E0 = (404.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,230.356,230.356,230.356,230.357,230.357],'cm^-1')),
        HinderedRotor(inertia=(0.547526,'amu*angstrom^2'), symmetry=1, barrier=(20.6174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547526,'amu*angstrom^2'), symmetry=1, barrier=(20.6173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547526,'amu*angstrom^2'), symmetry=1, barrier=(20.6173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205763,0.0705641,-4.53615e-05,-4.85405e-09,1.06225e-11,48825.4,31.4522], Tmin=(100,'K'), Tmax=(946.796,'K')), NASAPolynomial(coeffs=[21.3117,0.0127396,-3.40786e-06,5.70219e-10,-4.22333e-14,43423.9,-76.637], Tmin=(946.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O][CH]CC[C]=C=C[O](25713)',
    structure = SMILES('[O][CH]CC[C]=C=C[O]'),
    E0 = (516.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,356.816,357.226,358.22,358.995,360.402],'cm^-1')),
        HinderedRotor(inertia=(0.109973,'amu*angstrom^2'), symmetry=1, barrier=(9.94373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108458,'amu*angstrom^2'), symmetry=1, barrier=(9.92917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223508,'amu*angstrom^2'), symmetry=1, barrier=(20.544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273297,0.087318,-0.000119631,9.07142e-08,-2.78209e-11,62305.7,31.424], Tmin=(100,'K'), Tmax=(795.561,'K')), NASAPolynomial(coeffs=[11.3406,0.0316734,-1.47161e-05,2.7978e-09,-1.9398e-13,60544.7,-19.4381], Tmin=(795.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=COJ) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C[CH][CH]C=C[O](25714)',
    structure = SMILES('[O]C=C[CH][CH]C=C[O]'),
    E0 = (211.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,427.759,427.762,427.763,427.763,427.769],'cm^-1')),
        HinderedRotor(inertia=(0.241268,'amu*angstrom^2'), symmetry=1, barrier=(31.327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241248,'amu*angstrom^2'), symmetry=1, barrier=(31.3276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241253,'amu*angstrom^2'), symmetry=1, barrier=(31.3276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.670331,0.0508793,2.23081e-05,-8.07277e-08,3.87209e-11,25549.8,26.6246], Tmin=(100,'K'), Tmax=(937.52,'K')), NASAPolynomial(coeffs=[23.8372,0.00866731,-7.62591e-07,1.08997e-10,-1.63452e-14,18717.1,-96.921], Tmin=(937.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = '[O]C=C=[C]C[CH]C[O](25715)',
    structure = SMILES('[O]C=C=[C]C[CH]C[O]'),
    E0 = (536.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,180,180,180,663.043],'cm^-1')),
        HinderedRotor(inertia=(0.006901,'amu*angstrom^2'), symmetry=1, barrier=(2.15303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0432234,'amu*angstrom^2'), symmetry=1, barrier=(13.4843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0432227,'amu*angstrom^2'), symmetry=1, barrier=(13.4843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82801,0.0722838,-7.59232e-05,4.36149e-08,-1.03075e-11,64646.4,32.3809], Tmin=(100,'K'), Tmax=(1012.66,'K')), NASAPolynomial(coeffs=[11.4757,0.0302259,-1.36259e-05,2.60306e-09,-1.82888e-13,62489.9,-19.1222], Tmin=(1012.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=COJ) + radical(Cds_S) + radical(CCJCO)"""),
)

species(
    label = '[O]C=C[CH]C=[C]C[O](25716)',
    structure = SMILES('[O]C=C[CH]C=[C]C[O]'),
    E0 = (397.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,481.374,481.375,481.375,481.375,481.375],'cm^-1')),
        HinderedRotor(inertia=(0.141076,'amu*angstrom^2'), symmetry=1, barrier=(23.1979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00477162,'amu*angstrom^2'), symmetry=1, barrier=(42.6309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499889,'amu*angstrom^2'), symmetry=1, barrier=(82.1991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970496,0.0575368,-3.17641e-05,-6.68307e-10,4.59864e-12,47928.8,31.5327], Tmin=(100,'K'), Tmax=(1024.53,'K')), NASAPolynomial(coeffs=[14.0243,0.0243849,-9.30702e-06,1.68916e-09,-1.17642e-13,44319.1,-36.3233], Tmin=(1024.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJC=C) + radical(Cds_S) + radical(C=COJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C=[C]CC[CH][C]=O(25717)',
    structure = SMILES('[O]C=[C]CC[CH][C]=O'),
    E0 = (351.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1855,455,950,3025,407.5,1350,352.5,412.463,412.516,412.553,412.598],'cm^-1')),
        HinderedRotor(inertia=(0.000991109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1139,'amu*angstrom^2'), symmetry=1, barrier=(13.7534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113902,'amu*angstrom^2'), symmetry=1, barrier=(13.7539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17262,'amu*angstrom^2'), symmetry=1, barrier=(20.8318,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258598,0.0746736,-7.57166e-05,3.92101e-08,-7.9475e-12,42371.2,32.3507], Tmin=(100,'K'), Tmax=(1209.81,'K')), NASAPolynomial(coeffs=[17.1288,0.0188951,-6.55817e-06,1.10004e-09,-7.22071e-14,38289.2,-52.2514], Tmin=(1209.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCCJ=O) + radical(Cds_S) + radical(CCJCHO)"""),
)

species(
    label = '[O]C=[C]CC=[C]C[O](25718)',
    structure = SMILES('[O]C=[C]CC=[C]C[O]'),
    E0 = (533.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,369.492,369.492,369.492,369.492,2334.22],'cm^-1')),
        HinderedRotor(inertia=(0.134871,'amu*angstrom^2'), symmetry=1, barrier=(13.0665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134871,'amu*angstrom^2'), symmetry=1, barrier=(13.0665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134871,'amu*angstrom^2'), symmetry=1, barrier=(13.0665,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0484,0.069249,-7.47924e-05,4.69463e-08,-1.24427e-11,64285.5,32.0716], Tmin=(100,'K'), Tmax=(898.326,'K')), NASAPolynomial(coeffs=[9.01394,0.0337802,-1.55669e-05,2.99322e-09,-2.10626e-13,62854.4,-5.50327], Tmin=(898.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CCOJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C[CH]C[CH][C]=O(25719)',
    structure = SMILES('[O]C=C[CH]C[CH][C]=O'),
    E0 = (254.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3050,390,425,1340,1360,335,370,504.422,504.779,504.822,504.909],'cm^-1')),
        HinderedRotor(inertia=(0.000662447,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102259,'amu*angstrom^2'), symmetry=1, barrier=(18.476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102026,'amu*angstrom^2'), symmetry=1, barrier=(18.4699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.484803,'amu*angstrom^2'), symmetry=1, barrier=(87.6398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626352,0.0632617,-3.65639e-05,-5.30398e-09,8.7405e-12,30727.5,29.4495], Tmin=(100,'K'), Tmax=(961.382,'K')), NASAPolynomial(coeffs=[17.5708,0.0182187,-6.00516e-06,1.04875e-09,-7.39836e-14,26293,-57.7488], Tmin=(961.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCHO) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = '[O][CH]CC[CH][C]=C=O(25720)',
    structure = SMILES('[O][CH]CC[CH][C]=C=O'),
    E0 = (446.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2120,512.5,787.5,3000,3050,390,425,1340,1360,335,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333745,0.0923768,-0.000153906,1.43959e-07,-5.21338e-11,53824,32.0192], Tmin=(100,'K'), Tmax=(841.962,'K')), NASAPolynomial(coeffs=[5.67319,0.0424001,-2.10254e-05,4.02838e-09,-2.76864e-13,53797.2,12.3581], Tmin=(841.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCsJOH) + radical(CCJC(C)=C=O) + radical(CCOJ) + radical(CCCJ=C=O)"""),
)

species(
    label = '[O]C[C]=CC[CH][C]=O(25721)',
    structure = SMILES('[O]C[C]=CC[CH][C]=O'),
    E0 = (478.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1855,455,950,3025,407.5,1350,352.5,180,1352.72,1352.78,1352.81],'cm^-1')),
        HinderedRotor(inertia=(0.282739,'amu*angstrom^2'), symmetry=1, barrier=(6.50072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282675,'amu*angstrom^2'), symmetry=1, barrier=(6.49926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282735,'amu*angstrom^2'), symmetry=1, barrier=(6.50063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282715,'amu*angstrom^2'), symmetry=1, barrier=(6.50018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.977897,0.0752012,-0.000113134,1.03804e-07,-3.76584e-11,57683.6,32.6754], Tmin=(100,'K'), Tmax=(840.894,'K')), NASAPolynomial(coeffs=[4.44475,0.0405871,-1.90614e-05,3.59336e-09,-2.45691e-13,57741.3,20.3604], Tmin=(840.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(Cds_S) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = '[O]C[CH]C[CH][C]=C=O(25722)',
    structure = SMILES('[O]C[CH]C[CH][C]=C=O'),
    E0 = (466.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2120,512.5,787.5,3000,3050,390,425,1340,1360,335,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790991,0.0785032,-0.000114408,1.02812e-07,-3.75283e-11,56169,33.3258], Tmin=(100,'K'), Tmax=(801.835,'K')), NASAPolynomial(coeffs=[5.58185,0.0413714,-2.01917e-05,3.89697e-09,-2.71203e-13,55826.1,13.9232], Tmin=(801.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCOJ) + radical(CCCJ=C=O) + radical(CCJCO) + radical(CCJC(C)=C=O)"""),
)

species(
    label = 'C#CC([O])C([O])C=C(22582)',
    structure = SMILES('C#CC([O])C([O])C=C'),
    E0 = (344.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3010,987.5,1337.5,450,1655,408.426,408.431,408.455,408.486],'cm^-1')),
        HinderedRotor(inertia=(0.191236,'amu*angstrom^2'), symmetry=1, barrier=(22.6485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191292,'amu*angstrom^2'), symmetry=1, barrier=(22.6488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191389,'amu*angstrom^2'), symmetry=1, barrier=(22.6493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424544,0.0678809,-4.67261e-05,3.22755e-09,6.10177e-12,41603.2,31.6407], Tmin=(100,'K'), Tmax=(966.744,'K')), NASAPolynomial(coeffs=[18.6193,0.0169955,-5.62688e-06,9.90032e-10,-7.02363e-14,36945.2,-61.4194], Tmin=(966.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C=C1[CH]C[CH]C1[O](25723)',
    structure = SMILES('[O]C=C1[CH]C[CH]C1[O]'),
    E0 = (324.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5762,0.023065,9.95423e-05,-1.56573e-07,6.41213e-11,39128.2,26.1099], Tmin=(100,'K'), Tmax=(945.658,'K')), NASAPolynomial(coeffs=[22.8419,0.00985879,-1.242e-06,2.9577e-10,-3.69161e-14,31674.6,-93.44], Tmin=(945.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(methylenecyclopentane) + radical(C=COJ) + radical(Allyl_S) + radical(CC(C)OJ) + radical(CCJCO)"""),
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
    label = '[CH]=CCC=C=C[O](24301)',
    structure = SMILES('[CH]=CCC=C=C[O]'),
    E0 = (409.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10608,'amu*angstrom^2'), symmetry=1, barrier=(25.4309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1054,'amu*angstrom^2'), symmetry=1, barrier=(25.4154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916992,0.0557548,-2.4102e-05,-1.6588e-08,1.24936e-11,49324.6,26.4213], Tmin=(100,'K'), Tmax=(961.807,'K')), NASAPolynomial(coeffs=[17.6378,0.0147721,-4.72217e-06,8.48375e-10,-6.23169e-14,44787.3,-60.4619], Tmin=(961.807,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = 'C#C[CH]CC=C[O](20022)',
    structure = SMILES('C#C[CH]CC=C[O]'),
    E0 = (308.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,3025,407.5,1350,352.5,287.209,289.294,289.85],'cm^-1')),
        HinderedRotor(inertia=(0.47549,'amu*angstrom^2'), symmetry=1, barrier=(27.955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59061,'amu*angstrom^2'), symmetry=1, barrier=(34.3175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12682,'amu*angstrom^2'), symmetry=1, barrier=(66.5005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0647921,0.0677002,-6.37812e-05,3.0016e-08,-5.2724e-12,37291.1,28.5977], Tmin=(100,'K'), Tmax=(1651.61,'K')), NASAPolynomial(coeffs=[17.2967,0.0123796,-1.48396e-06,4.11297e-12,6.9518e-15,33366.6,-58.3917], Tmin=(1651.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=COJ)"""),
)

species(
    label = '[O]C=[C]C1CC1C=O(25368)',
    structure = SMILES('[O]C=[C]C1CC1C=O'),
    E0 = (172.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976809,0.0501206,4.52678e-06,-4.98626e-08,2.48463e-11,20870,27.67], Tmin=(100,'K'), Tmax=(953.482,'K')), NASAPolynomial(coeffs=[18.7955,0.0156498,-4.61388e-06,8.35914e-10,-6.39111e-14,15641,-67.0478], Tmin=(953.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C=CC=CC=O(25724)',
    structure = SMILES('[O]C=C=CC=CC=O'),
    E0 = (43.6819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20709,'amu*angstrom^2'), symmetry=1, barrier=(27.7535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20965,'amu*angstrom^2'), symmetry=1, barrier=(27.8123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347232,0.0691822,-5.53387e-05,1.42945e-08,1.13954e-12,5394.91,26.3143], Tmin=(100,'K'), Tmax=(1045.7,'K')), NASAPolynomial(coeffs=[19.6858,0.0151103,-6.32322e-06,1.24573e-09,-9.19903e-14,262.33,-73.0499], Tmin=(1045.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.6819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[O]C=[C]C=CCC=O(25406)',
    structure = SMILES('[O]C=[C]C=CCC=O'),
    E0 = (90.2986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547981,0.0613444,-2.36644e-05,-1.98172e-08,1.32361e-11,10997.6,27.5579], Tmin=(100,'K'), Tmax=(1001.62,'K')), NASAPolynomial(coeffs=[19.5092,0.0179589,-7.11818e-06,1.40209e-09,-1.05183e-13,5577.12,-72.0471], Tmin=(1001.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.2986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CCC[C]=O(22577)',
    structure = SMILES('[O]C=C=CCC[C]=O'),
    E0 = (110.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.727067,'amu*angstrom^2'), symmetry=1, barrier=(16.7167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.727233,'amu*angstrom^2'), symmetry=1, barrier=(16.7205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.727536,'amu*angstrom^2'), symmetry=1, barrier=(16.7275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.388541,0.074216,-7.46134e-05,3.83428e-08,-7.77991e-12,13377.4,30.2659], Tmin=(100,'K'), Tmax=(1199.88,'K')), NASAPolynomial(coeffs=[16.3039,0.0211601,-8.2876e-06,1.49189e-09,-1.01943e-13,9558.03,-49.4167], Tmin=(1199.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCCJ=O) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=[C]CCC=O(25725)',
    structure = SMILES('[O]C=C=[C]CCC=O'),
    E0 = (187.989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,181.407,181.407,181.407],'cm^-1')),
        HinderedRotor(inertia=(0.717664,'amu*angstrom^2'), symmetry=1, barrier=(16.7593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717664,'amu*angstrom^2'), symmetry=1, barrier=(16.7593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717666,'amu*angstrom^2'), symmetry=1, barrier=(16.7593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684673,0.0727976,-7.35687e-05,3.92504e-08,-8.47662e-12,22729.4,28.861], Tmin=(100,'K'), Tmax=(1112.29,'K')), NASAPolynomial(coeffs=[13.4203,0.026998,-1.18045e-05,2.23112e-09,-1.56098e-13,19896.3,-33.9362], Tmin=(1112.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'O=C=[C][CH]CCC=O(25726)',
    structure = SMILES('O=C=[C][CH]CCC=O'),
    E0 = (117.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,3025,407.5,1350,352.5,207.133,210.76,2491.86],'cm^-1')),
        HinderedRotor(inertia=(0.358097,'amu*angstrom^2'), symmetry=1, barrier=(10.8636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00386816,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30532,'amu*angstrom^2'), symmetry=1, barrier=(40.2117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32947,'amu*angstrom^2'), symmetry=1, barrier=(40.2245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54823,0.0648669,-3.85194e-05,-5.53219e-08,7.63729e-11,14213.8,26.7624], Tmin=(100,'K'), Tmax=(474.342,'K')), NASAPolynomial(coeffs=[6.88885,0.0392685,-1.90372e-05,3.68581e-09,-2.57956e-13,13488.4,2.67514], Tmin=(474.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CCCJ=C=O) + radical(CCJC(C)=C=O)"""),
)

species(
    label = '[O]C=C1[CH]CC1C=O(25727)',
    structure = SMILES('[O]C=C1[CH]CC1C=O'),
    E0 = (60.6614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.6648,-0.0113478,0.000123746,-1.21222e-07,2.80117e-11,6937.54,-17.8616], Tmin=(100,'K'), Tmax=(1722.81,'K')), NASAPolynomial(coeffs=[77.5207,0.0219503,-6.93877e-05,1.70312e-08,-1.26746e-12,-44076,-457.967], Tmin=(1722.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.6614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(methylenecyclobutane) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = '[O]C1[C]=CCC1C=O(25478)',
    structure = SMILES('[O]C1[C]=CCC1C=O'),
    E0 = (224.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39804,0.0438242,6.78485e-06,-3.90307e-08,1.73519e-11,27132.7,26.3071], Tmin=(100,'K'), Tmax=(1014.2,'K')), NASAPolynomial(coeffs=[14.2224,0.024136,-9.78381e-06,1.89201e-09,-1.38269e-13,22942.7,-43.5762], Tmin=(1014.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(cyclopentene-vinyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'C1=CC[CH][CH]OOC=1(25593)',
    structure = SMILES('C1=CC[CH][CH]OOC=1'),
    E0 = (508.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80419,0.0304462,5.06067e-05,-8.29643e-08,3.20905e-11,61261.5,23.9613], Tmin=(100,'K'), Tmax=(1000.54,'K')), NASAPolynomial(coeffs=[13.4392,0.0278692,-1.14014e-05,2.24278e-09,-1.66546e-13,56733.9,-43.1678], Tmin=(1000.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclooctane) + radical(CCJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = 'O=CC=CC=C=CO(25728)',
    structure = SMILES('O=CC=CC=C=CO'),
    E0 = (-97.7808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0434416,0.0759566,-5.83246e-05,7.98723e-09,5.69611e-12,-11603.1,26.3036], Tmin=(100,'K'), Tmax=(979.023,'K')), NASAPolynomial(coeffs=[22.5942,0.0119884,-4.01707e-06,7.64445e-10,-5.83952e-14,-17402.6,-89.4112], Tmin=(979.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.7808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O=C=C=CCCC=O(25729)',
    structure = SMILES('O=C=C=CCCC=O'),
    E0 = (-5.75617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0763,0.0438425,-2.8053e-05,8.947e-09,-1.18728e-12,-625.082,10.5808], Tmin=(100,'K'), Tmax=(1641.85,'K')), NASAPolynomial(coeffs=[9.65986,0.0253671,-1.11739e-05,2.09339e-09,-1.43708e-13,-3115.32,-29.7657], Tmin=(1641.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.75617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C(C=O)C=C=C[O](22596)',
    structure = SMILES('[CH2]C(C=O)C=C=C[O]'),
    E0 = (153.763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.858345,'amu*angstrom^2'), symmetry=1, barrier=(19.735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85823,'amu*angstrom^2'), symmetry=1, barrier=(19.7324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858549,'amu*angstrom^2'), symmetry=1, barrier=(19.7397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4311.52,'J/mol'), sigma=(6.83844,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=673.45 K, Pc=30.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.362348,0.0800131,-8.91188e-05,5.1578e-08,-1.18964e-11,18624.5,29.091], Tmin=(100,'K'), Tmax=(1053.85,'K')), NASAPolynomial(coeffs=[14.8601,0.0249841,-1.07915e-05,2.02689e-09,-1.41383e-13,15568.8,-41.6122], Tmin=(1053.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CJC(C)C=O)"""),
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
    label = '[CH]CC=C=C[O](22690)',
    structure = SMILES('[CH]CC=C=C[O]'),
    E0 = (504.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,440.804,441.131,441.44,441.765,443.303],'cm^-1')),
        HinderedRotor(inertia=(0.139422,'amu*angstrom^2'), symmetry=1, barrier=(19.2396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140211,'amu*angstrom^2'), symmetry=1, barrier=(19.2388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35618,0.0474549,-1.98078e-05,-1.658e-08,1.20553e-11,60835.4,22.7003], Tmin=(100,'K'), Tmax=(950.147,'K')), NASAPolynomial(coeffs=[16.2492,0.0108067,-3.07525e-06,5.34262e-10,-3.98986e-14,56829.5,-54.5757], Tmin=(950.147,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C=CCC#CC=O(25730)',
    structure = SMILES('[O]C=CCC#CC=O'),
    E0 = (60.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2100,2250,500,550,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.09816,'amu*angstrom^2'), symmetry=1, barrier=(25.2489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09887,'amu*angstrom^2'), symmetry=1, barrier=(25.2652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09872,'amu*angstrom^2'), symmetry=1, barrier=(25.2618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.854715,0.0588036,-3.55188e-05,-8.83957e-10,5.66653e-12,7439.6,27.3767], Tmin=(100,'K'), Tmax=(1008.82,'K')), NASAPolynomial(coeffs=[16.5072,0.0181707,-6.96587e-06,1.30389e-09,-9.38054e-14,3191.03,-53.6797], Tmin=(1008.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CC[C]=CC=O(25731)',
    structure = SMILES('[O]C=CC[C]=CC=O'),
    E0 = (135.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.795126,'amu*angstrom^2'), symmetry=1, barrier=(18.2815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.795313,'amu*angstrom^2'), symmetry=1, barrier=(18.2858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794299,'amu*angstrom^2'), symmetry=1, barrier=(18.2625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582611,0.068327,-5.86277e-05,2.45842e-08,-4.08711e-12,16458.3,29.9183], Tmin=(100,'K'), Tmax=(1436.45,'K')), NASAPolynomial(coeffs=[17.4345,0.0214009,-9.62585e-06,1.84222e-09,-1.29125e-13,11616.9,-57.4857], Tmin=(1436.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=CC[CH]C=C=O(25732)',
    structure = SMILES('[O]C=CC[CH]C=C=O'),
    E0 = (58.5105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440949,0.0691255,-5.74624e-05,2.10545e-08,-2.14756e-12,7173.13,30.5599], Tmin=(100,'K'), Tmax=(1087.13,'K')), NASAPolynomial(coeffs=[17.001,0.0208028,-8.18467e-06,1.50413e-09,-1.05058e-13,2827.49,-54.1428], Tmin=(1087.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.5105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCJC(C)=C=O) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C[CH]C=CC=O(25733)',
    structure = SMILES('[O]C=C[CH]C=CC=O'),
    E0 = (-0.343332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885415,0.0520678,3.76489e-07,-4.34835e-08,2.14359e-11,85.5156,27.3267], Tmin=(100,'K'), Tmax=(984.157,'K')), NASAPolynomial(coeffs=[18.7787,0.0180537,-6.78241e-06,1.33347e-09,-1.01403e-13,-5311.16,-68.2368], Tmin=(984.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.343332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = '[O]C=[C]CC=CC=O(25734)',
    structure = SMILES('[O]C=[C]CC=CC=O'),
    E0 = (135.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.795126,'amu*angstrom^2'), symmetry=1, barrier=(18.2815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.795313,'amu*angstrom^2'), symmetry=1, barrier=(18.2858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794299,'amu*angstrom^2'), symmetry=1, barrier=(18.2625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582611,0.068327,-5.86277e-05,2.45842e-08,-4.08711e-12,16458.3,29.9183], Tmin=(100,'K'), Tmax=(1436.45,'K')), NASAPolynomial(coeffs=[17.4345,0.0214009,-9.62585e-06,1.84222e-09,-1.29125e-13,11616.9,-57.4857], Tmin=(1436.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'O=[C][CH]CC=CC=O(25735)',
    structure = SMILES('O=[C][CH]CC=CC=O'),
    E0 = (80.9006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31604,0.0644673,-6.11839e-05,3.30175e-08,-7.72974e-12,9822.21,27.6616], Tmin=(100,'K'), Tmax=(988.605,'K')), NASAPolynomial(coeffs=[8.46561,0.0355392,-1.72913e-05,3.4183e-09,-2.44582e-13,8408.6,-6.74888], Tmin=(988.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.9006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = '[O]C1[CH]CC=C1C=O(25736)',
    structure = SMILES('[O]C1[CH]CC=C1C=O'),
    E0 = (148.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52074,0.0371962,3.03895e-05,-6.45157e-08,2.59944e-11,17982.6,27.0711], Tmin=(100,'K'), Tmax=(1019.64,'K')), NASAPolynomial(coeffs=[15.7066,0.023183,-1.02492e-05,2.10457e-09,-1.59509e-13,12925.2,-52.2573], Tmin=(1019.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopentene) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C=CCC=CC=O(25737)',
    structure = SMILES('O=C=CCC=CC=O'),
    E0 = (-121.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06311,0.0713225,-8.56804e-05,6.4027e-08,-2.08636e-11,-14480,25.7232], Tmin=(100,'K'), Tmax=(729.242,'K')), NASAPolynomial(coeffs=[6.8815,0.0394071,-2.00311e-05,4.00969e-09,-2.87902e-13,-15328.6,-0.509874], Tmin=(729.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'O=CC1=CCC=CO1(25521)',
    structure = SMILES('O=CC1=CCC=CO1'),
    E0 = (-146.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13876,0.0470674,8.11651e-06,-4.53035e-08,2.02806e-11,-17487.4,23.0265], Tmin=(100,'K'), Tmax=(1014.62,'K')), NASAPolynomial(coeffs=[16.7844,0.0221371,-9.35811e-06,1.87723e-09,-1.40748e-13,-22553.9,-62.0034], Tmin=(1014.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(1,4-Cyclohexadiene)"""),
)

species(
    label = '[C]=CCC=C[O](13810)',
    structure = SMILES('[C]=CCC=C[O]'),
    E0 = (579.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,193.415,193.417,193.421],'cm^-1')),
        HinderedRotor(inertia=(0.754641,'amu*angstrom^2'), symmetry=1, barrier=(20.0338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754626,'amu*angstrom^2'), symmetry=1, barrier=(20.0339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56853,0.0433028,-1.47111e-05,-1.78279e-08,1.14878e-11,69797.1,23.1981], Tmin=(100,'K'), Tmax=(964.286,'K')), NASAPolynomial(coeffs=[14.8673,0.0122777,-4.00083e-06,7.28701e-10,-5.3868e-14,66110,-46.2965], Tmin=(964.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(579.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'O=CC=CC=CC=O(25738)',
    structure = SMILES('O=CC=CC=CC=O'),
    E0 = (-153.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44059,0.062269,-4.82631e-05,1.79372e-08,-2.76351e-12,-18322,22.6304], Tmin=(100,'K'), Tmax=(1443.7,'K')), NASAPolynomial(coeffs=[12.3933,0.0319227,-1.67332e-05,3.37736e-09,-2.42232e-13,-21484.4,-34.2319], Tmin=(1443.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'O=CC1=CCC1C=O(25363)',
    structure = SMILES('O=CC1=CCC1C=O'),
    E0 = (-99.8132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24895,0.0537625,-3.18406e-05,7.78255e-09,-5.84595e-13,-11900.1,24.6169], Tmin=(100,'K'), Tmax=(1630.24,'K')), NASAPolynomial(coeffs=[16.6641,0.0234639,-1.08857e-05,2.04445e-09,-1.3881e-13,-17926,-60.3525], Tmin=(1630.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.8132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene)"""),
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
    E0 = (94.6684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (174.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (249.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (290.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (353.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (355.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (344.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (341.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (209.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (490.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (224.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (156.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (342.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (201.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (218.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (447.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (361.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (408.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (544.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (544.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (489.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (472.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (308.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (281.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (220.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (207.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (237.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (276.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (119.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (119.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (435.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (335.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (399.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (332.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (497.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (605.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (300.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (575.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (436.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (368.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (572.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (279.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (471.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (503.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (491.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (474.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (815.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (715.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (202.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (264.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (193.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (201.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (266.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (329.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (180.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (219.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (237.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (508.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (112.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (133.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (313.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (572.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (282.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (216.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (332.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (257.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (241.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (530.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (127.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (148.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (133.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (102.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (647.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (187.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (102.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['C=CC=O(5269)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=[C]C1CC=CO1(25557)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.3435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=CC1C[CH][C]=CO1(25456)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.04811e+11,'s^-1'), n=0.222, Ea=(155.084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra;radadd_intra_O] + [R7_SMMS;doublebond_intra;radadd_intra] for rate rule [R7_SMMS;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[O]C=C=CCC=C=O(25685)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[O]C=CCC=C=C=O(25686)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C=C=CCC=[C]O(25687)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=CC[CH]C#CO(25688)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=C=C[CH]C=CO(25690)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.92799e+06,'s^-1'), n=1.7075, Ea=(114.432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C=CCC#C[CH]O(25691)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C=C=[C]CC=CO(25692)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSMS;Cd_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=C[CH]C=C=CO(25693)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(60205.5,'s^-1'), n=1.86417, Ea=(61.9987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_SMMS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C=[C]CC=C=CO(25694)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R6H;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C=[C][CH]CC=CO(25695)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.21847e+06,'s^-1'), n=1.22418, Ea=(82.4275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;XH_out] for rate rule [R7H;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C][CH]CC=C=CO(25696)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.21847e+06,'s^-1'), n=1.22418, Ea=(82.4275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;XH_out] for rate rule [R7H;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C[O](602)', 'C=C[C]=C[O](22438)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.169e+08,'m^3/(mol*s)'), n=-0.455312, Ea=(0.377199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_pri_rad;Cd_rad] + [C_rad/H2/Cd;Y_rad] for rate rule [C_rad/H2/Cd;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C=C[O](5266)', '[CH]=C=C[O](8556)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_allenic;C_pri_rad] for rate rule [Cd_allenic;C_rad/H2/Cd]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[O]C=C=C[CH]C=C[O](25697)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.14e-13,'cm^3/(molecule*s)'), n=0.611, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 53 used for C_rad/H/CdCd;H_rad
Exact match found for rate rule [C_rad/H/CdCd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.8 to -1.8 kJ/mol.
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[O]C=[C]CC=C=C[O](25698)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[O]C=C=[C]CC=C[O](25699)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[O]C=C=CC[CH][C]=O(25700)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[O]C=CC[CH][C]=C=O(25701)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=C=CCC1[CH]O1(25702)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=CCC=C1[CH]O1(25703)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=CCC1[C]=CO1(25704)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=C1[CH]CC=CO1(25504)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.68243e+09,'s^-1'), n=0.4695, Ea=(113.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C1[C]=CCC=CO1(25636)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.22e+12,'s^-1'), n=-0.622, Ea=(142.884,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_CdCdd;radadd_intra] for rate rule [R7;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C1[CH]CC=C=CO1(25705)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.51044e+06,'s^-1'), n=0.97919, Ea=(181.92,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 179.3 to 181.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=C=CCC=C=CO(25706)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=C=C=CCC=CO(25707)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O][CH]CC=C[C]=C[O](25708)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C=[C]C[CH]C=C[O](25709)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C=C[C]CC=C[O](25710)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O][C]=C[CH]CC=C[O](25711)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C=[C]CC[C]=C[O](25712)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.16e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O][CH]CC[C]=C=C[O](25713)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C=C[CH][CH]C=C[O](25714)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C=C=[C]C[CH]C[O](25715)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]C=C[CH]C=[C]C[O](25716)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]C=[C]CC[CH][C]=O(25717)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O]C=[C]CC=[C]C[O](25718)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]C=C[CH]C[CH][C]=O(25719)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O][CH]CC[CH][C]=C=O(25720)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[O]C[C]=CC[CH][C]=O(25721)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O]C[CH]C[CH][C]=C=O(25722)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C#CC([O])C([O])C=C(22582)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.214e+09,'s^-1'), n=0.749, Ea=(200.242,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2 used for hex_1_ene_5_yne
Exact match found for rate rule [hex_1_ene_5_yne]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[O]C=C1[CH]C[CH]C1[O](25723)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction48',
    reactants = ['O(T)(63)', '[CH]=CCC=C=C[O](24301)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['O(T)(63)', 'C#C[CH]CC=C[O](20022)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=[C]C1CC1C=O(25368)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHDe]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['H(8)', '[O]C=C=CC=CC=O(25724)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(30.9498,'m^3/(mol*s)'), n=1.77449, Ea=(9.26302,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-OneDeH;HJ] for rate rule [Cds-CdH_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=CC=O(5269)', '[CH]=C=C[O](8556)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(0.0102751,'m^3/(mol*s)'), n=2.40501, Ea=(4.48561,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDeH;CJ] for rate rule [Cds-HH_Cds-COH;CdsJ=Cdd]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=[C]C=CCC=O(25406)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(8.03429e+08,'s^-1'), n=1.36289, Ea=(107.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[O]C=C=CCC[C]=O(22577)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[O]C=C=[C]CCC=O(25725)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_double;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_double;Cs_H_out_H/CO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['O=C=[C][CH]CCC=O(25726)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(59826.2,'s^-1'), n=1.86494, Ea=(62.4776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5H;Cd_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=C1[CH]CC1C=O(25727)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(3.78932e+07,'s^-1'), n=1.19089, Ea=(125.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHCO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C1[C]=CCC1C=O(25478)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.77218e+12,'s^-1'), n=-0.0488474, Ea=(143.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5;doublebond_intra;radadd_intra_csHCO] + [R5_linear;doublebond_intra;radadd_intra_csHDe] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_csHCO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['C1=CC[CH][CH]OOC=1(25593)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(413.9,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 409.5 to 413.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction60',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=CC=CC=C=CO(25728)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=C=C=CCCC=O(25729)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH]=O(373)', '[CH]CC=C=C[O](22690)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction64',
    reactants = ['H(8)', '[O]C=CCC#CC=O(25730)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]C=C[O](5266)', 'C#CC=O(21959)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(0.0835733,'m^3/(mol*s)'), n=2.41, Ea=(41.6726,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CsJ-CdHH] for rate rule [Ct-H_Ct-CO;CsJ-CdHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[O]C=CC[C]=CC=O(25731)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=CC[CH]C=C=O(25732)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C=C[CH]C=CC=O(25733)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe;Cs_H_out_H/Cd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[O]C=[C]CC=CC=O(25734)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=[C][CH]CC=CC=O(25735)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['[O]C1[CH]CC=C1C=O(25736)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(1.80657e+08,'s^-1'), n=0.835, Ea=(53.979,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 50.3 to 54.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction72',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=C=CCC=CC=O(25737)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=CC1=CCC=CO1(25521)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;CdsingleDe_rad_out;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction74',
    reactants = ['[CH]=O(373)', '[C]=CCC=C[O](13810)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=CC=CC=CC=O(25738)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[O]C=C=CCC=C[O](22588)'],
    products = ['O=CC1=CCC1C=O(25363)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriDe_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

network(
    label = '4889',
    isomers = [
        '[O]C=C=CCC=C[O](22588)',
    ],
    reactants = [
        ('C=CC=O(5269)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4889',
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

