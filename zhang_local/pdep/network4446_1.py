species(
    label = '[CH][CH]C=C=[CH](19252)',
    structure = SMILES('[CH][CH]C=C=[CH]'),
    E0 = (867.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,524.185,563.719,564.949],'cm^-1')),
        HinderedRotor(inertia=(0.0200302,'amu*angstrom^2'), symmetry=1, barrier=(79.5804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354495,'amu*angstrom^2'), symmetry=1, barrier=(79.5766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19239,0.0344402,-1.76598e-05,-4.40435e-09,5.18614e-12,104450,17.9495], Tmin=(100,'K'), Tmax=(942.255,'K')), NASAPolynomial(coeffs=[10.2948,0.0130216,-4.22214e-06,7.05005e-10,-4.76148e-14,102347,-23.7152], Tmin=(942.255,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(867.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CCJ2_triplet) + radical(Allyl_S)"""),
)

species(
    label = 'C#C(582)',
    structure = SMILES('C#C'),
    E0 = (214.792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,559.488,618.58,3890.62],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03575,0.00771239,2.53493e-06,-1.08133e-08,5.50757e-12,25852.6,4.54461], Tmin=(100,'K'), Tmax=(888.627,'K')), NASAPolynomial(coeffs=[5.76205,0.00237159,-1.49583e-07,-2.19155e-11,2.21779e-15,25094.5,-9.82608], Tmin=(888.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.792,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Ct-CtH) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C=[CH](18734)',
    structure = SMILES('[CH]=C=[CH]'),
    E0 = (491.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,239.877,511.233,1743.98,1746.51,1747.6,1753.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.766,0.0170203,-1.57568e-05,7.95984e-09,-1.4265e-12,59188.9,11.2142], Tmin=(100,'K'), Tmax=(1806.04,'K')), NASAPolynomial(coeffs=[4.81405,0.00509933,2.77647e-07,-2.23082e-10,1.96202e-14,59653.5,3.45727], Tmin=(1806.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
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
    label = '[CH]C=C=C=[CH](21421)',
    structure = SMILES('[CH]C=C=C=[CH]'),
    E0 = (812.442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,563.333,586.667,610,1970,2140,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.10357,'amu*angstrom^2'), symmetry=1, barrier=(48.3651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07136,0.0400756,-3.45202e-05,1.69333e-08,-3.45179e-12,97785.6,17.5421], Tmin=(100,'K'), Tmax=(1167.18,'K')), NASAPolynomial(coeffs=[8.33549,0.0186082,-6.93152e-06,1.17533e-09,-7.65903e-14,96323.4,-13.6471], Tmin=(1167.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(812.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C[C]=C=[CH](21422)',
    structure = SMILES('[CH]C[C]=C=[CH]'),
    E0 = (964.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,1685,370,356.801,356.806,2406.53],'cm^-1')),
        HinderedRotor(inertia=(0.122399,'amu*angstrom^2'), symmetry=1, barrier=(11.0579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.79857,'amu*angstrom^2'), symmetry=1, barrier=(72.1419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99373,0.0438701,-4.9978e-05,3.14598e-08,-7.9315e-12,116086,20.2437], Tmin=(100,'K'), Tmax=(968.996,'K')), NASAPolynomial(coeffs=[8.93134,0.0152319,-5.64634e-06,9.59807e-10,-6.25778e-14,114741,-13.0078], Tmin=(968.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(964.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C=[C][C]=C(19273)',
    structure = SMILES('[CH]C=[C][C]=C'),
    E0 = (826.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15488,'amu*angstrom^2'), symmetry=1, barrier=(49.545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15199,'amu*angstrom^2'), symmetry=1, barrier=(49.4784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97318,0.0442752,-3.59148e-05,1.31176e-08,-1.0745e-13,99487,17.3868], Tmin=(100,'K'), Tmax=(793.849,'K')), NASAPolynomial(coeffs=[7.85963,0.0220176,-7.84619e-06,1.29259e-09,-8.28058e-14,98319.1,-11.1221], Tmin=(793.849,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(826.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=[C][C]=[CH](21423)',
    structure = SMILES('[CH]C=[C][C]=[CH]'),
    E0 = (1073.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13504,'amu*angstrom^2'), symmetry=1, barrier=(49.0888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13552,'amu*angstrom^2'), symmetry=1, barrier=(49.0999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75888,0.0496114,-5.97279e-05,4.11383e-08,-1.12717e-11,129213,18.5492], Tmin=(100,'K'), Tmax=(991.755,'K')), NASAPolynomial(coeffs=[8.35249,0.0187317,-6.54095e-06,1.02802e-09,-6.23706e-14,128116,-12.1439], Tmin=(991.755,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1073.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1[CH]C1=[CH](21424)',
    structure = SMILES('[CH]C1[CH]C1=[CH]'),
    E0 = (991.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45343,0.0211814,2.94054e-05,-5.84283e-08,2.52969e-11,119331,15.334], Tmin=(100,'K'), Tmax=(941.461,'K')), NASAPolynomial(coeffs=[13.1099,0.00835646,-1.86492e-06,3.274e-10,-2.76576e-14,115887,-43.0726], Tmin=(941.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(991.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][CH]C1[C]=C1(21425)',
    structure = SMILES('[CH][CH]C1[C]=C1'),
    E0 = (1088.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26514,0.0411483,-4.97357e-05,3.72328e-08,-1.18898e-11,130949,18.1593], Tmin=(100,'K'), Tmax=(751.151,'K')), NASAPolynomial(coeffs=[6.00782,0.0212181,-9.93684e-06,1.91056e-09,-1.33847e-13,130387,1.17395], Tmin=(751.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1088.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Cs_S) + radical(CCJ2_triplet) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = '[CH]=CC#C[CH2](19271)',
    structure = SMILES('[CH]=CC#C[CH2]'),
    E0 = (632.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2100,2250,500,550,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(1.49942,'amu*angstrom^2'), symmetry=1, barrier=(34.4747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49786,'amu*angstrom^2'), symmetry=1, barrier=(34.4387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37102,0.0279414,-1.39095e-07,-2.08295e-08,1.00915e-11,76081.2,17.6168], Tmin=(100,'K'), Tmax=(992.381,'K')), NASAPolynomial(coeffs=[10.6694,0.0127835,-4.87417e-06,9.24016e-10,-6.74612e-14,73533.5,-26.8926], Tmin=(992.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl) + radical(Cds_P)"""),
)

species(
    label = '[CH][CH](5289)',
    structure = SMILES('[CH][CH]'),
    E0 = (799.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([387.411,1315.31,1315.42,1316.38,1321.45],'cm^-1')),
        HinderedRotor(inertia=(0.00939611,'amu*angstrom^2'), symmetry=1, barrier=(0.481197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.68559,0.00274596,1.66519e-05,-2.51103e-08,1.03855e-11,96195.3,7.55883], Tmin=(100,'K'), Tmax=(914.567,'K')), NASAPolynomial(coeffs=[6.28769,0.00215456,-7.39604e-08,-1.89738e-11,6.135e-16,95268.1,-7.22932], Tmin=(914.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(799.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]=C[C]=[CH](21426)',
    structure = SMILES('[CH][C]=C[C]=[CH]'),
    E0 = (1112.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13683,'amu*angstrom^2'), symmetry=1, barrier=(49.1299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1316,'amu*angstrom^2'), symmetry=1, barrier=(49.0098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95039,0.0458062,-4.82344e-05,2.71146e-08,-5.46777e-12,133878,18.0544], Tmin=(100,'K'), Tmax=(785.352,'K')), NASAPolynomial(coeffs=[8.34208,0.0187805,-7.17602e-06,1.22516e-09,-7.99173e-14,132703,-12.3229], Tmin=(785.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1112.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=[C]C=C[CH](21427)',
    structure = SMILES('[C]=[C]C=C[CH]'),
    E0 = (1185.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.1391,'amu*angstrom^2'), symmetry=1, barrier=(49.1822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15895,'amu*angstrom^2'), symmetry=1, barrier=(49.6385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97287,0.0447068,-4.66913e-05,2.81817e-08,-7.02024e-12,142677,18.1628], Tmin=(100,'K'), Tmax=(967.746,'K')), NASAPolynomial(coeffs=[8.05196,0.0195802,-7.74572e-06,1.35291e-09,-8.95638e-14,141500,-10.9661], Tmin=(967.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1185.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CdCdJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=C=C[CH](21428)',
    structure = SMILES('[CH]C=C=C[CH]'),
    E0 = (852.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92581,0.040213,-1.89945e-05,2.83431e-09,2.21484e-13,102557,19.5094], Tmin=(100,'K'), Tmax=(1422.23,'K')), NASAPolynomial(coeffs=[8.89175,0.0265137,-1.07608e-05,1.88787e-09,-1.24216e-13,99979.2,-18.6459], Tmin=(1422.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(852.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C=[C][CH2](19274)',
    structure = SMILES('[CH]=[C]C=[C][CH2]'),
    E0 = (859.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1670,1700,300,440,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(2.0544,'amu*angstrom^2'), symmetry=1, barrier=(47.2346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05464,'amu*angstrom^2'), symmetry=1, barrier=(47.2402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93557,0.0405981,-3.87335e-05,2.02214e-08,-4.17345e-12,103499,19.2405], Tmin=(100,'K'), Tmax=(1265.13,'K')), NASAPolynomial(coeffs=[9.93801,0.0132358,-3.8481e-06,5.50841e-10,-3.19634e-14,101639,-20.5969], Tmin=(1265.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(859.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1[CH][C]=C1(21429)',
    structure = SMILES('[CH]C1[CH][C]=C1'),
    E0 = (977.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84395,0.0113557,5.19156e-05,-7.7351e-08,3.08154e-11,117676,14.4337], Tmin=(100,'K'), Tmax=(950.468,'K')), NASAPolynomial(coeffs=[11.499,0.0105363,-2.98218e-06,5.67601e-10,-4.60693e-14,114422,-35.3426], Tmin=(950.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(977.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CCJ2_triplet) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH]C=[CH](19276)',
    structure = SMILES('[CH]C=[CH]'),
    E0 = (623.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.17043,'amu*angstrom^2'), symmetry=1, barrier=(49.9024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24639,0.0115967,1.7963e-05,-2.77835e-08,1.05348e-11,75069.9,10.6897], Tmin=(100,'K'), Tmax=(978.982,'K')), NASAPolynomial(coeffs=[5.45376,0.0147826,-5.61896e-06,1.00989e-09,-7.00625e-14,74052.9,-2.8998], Tmin=(978.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=[CH](18830)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([562.459,1392.74,3112.83],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86574,0.0023643,1.31489e-06,-2.64796e-09,1.00932e-12,101918,6.62295], Tmin=(100,'K'), Tmax=(1016.96,'K')), NASAPolynomial(coeffs=[4.17915,0.00251347,-9.43473e-07,1.6872e-10,-1.15831e-14,101782,4.75427], Tmin=(1016.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C=C=[CH](21430)',
    structure = SMILES('[CH]=[C]C=C=[CH]'),
    E0 = (835.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,1685,370,540,610,2055,197.988,198.111,199.477,200.18],'cm^-1')),
        HinderedRotor(inertia=(0.0653758,'amu*angstrom^2'), symmetry=1, barrier=(71.8386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98287,0.0425038,-4.25398e-05,1.42591e-08,1.96062e-12,100510,16.0374], Tmin=(100,'K'), Tmax=(798.852,'K')), NASAPolynomial(coeffs=[11.6648,0.0071383,-7.57417e-07,-6.06055e-11,1.11936e-14,98544.8,-31.1168], Tmin=(798.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(835.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[CH](583)',
    structure = SMILES('[CH]=[CH]'),
    E0 = (536.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([637.691,1081.65,1081.98,1082.08,3058.36,3477.84],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.83395,-0.000554326,2.20867e-05,-2.90276e-08,1.14365e-11,64516.9,6.06922], Tmin=(100,'K'), Tmax=(916.167,'K')), NASAPolynomial(coeffs=[5.69903,0.00213261,-4.3877e-08,-2.1344e-11,5.56752e-16,63720.6,-5.24595], Tmin=(916.167,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]=[CH](21256)',
    structure = SMILES('[CH][C]=[CH]'),
    E0 = (861.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.1891,'amu*angstrom^2'), symmetry=1, barrier=(50.3317,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18317,0.0164338,-7.13252e-06,1.19383e-09,-3.27944e-14,103675,12.0918], Tmin=(100,'K'), Tmax=(1799.19,'K')), NASAPolynomial(coeffs=[6.32962,0.0112581,-4.33439e-06,7.19107e-10,-4.49321e-14,102248,-5.75439], Tmin=(1799.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(861.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C][C]=C[CH2](21431)',
    structure = SMILES('[CH]=[C][C]=C[CH2]'),
    E0 = (821.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14987,0.0393004,-3.07102e-05,6.09782e-09,3.53696e-12,98816.3,17.6059], Tmin=(100,'K'), Tmax=(793.488,'K')), NASAPolynomial(coeffs=[9.16594,0.0145565,-4.01894e-06,5.46848e-10,-3.06094e-14,97368.4,-16.7273], Tmin=(793.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(821.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1[CH]C=C1(21418)',
    structure = SMILES('[CH]=C1[CH]C=C1'),
    E0 = (568.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89236,0.000530591,0.000102751,-1.45866e-07,5.96435e-11,68442.3,11.8698], Tmin=(100,'K'), Tmax=(916.18,'K')), NASAPolynomial(coeffs=[17.3812,-0.00108824,4.48409e-06,-9.2893e-10,5.628e-14,63200.5,-70.88], Tmin=(916.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(C=CCJC=C)"""),
)

species(
    label = '[C][CH]C=C=[CH](21432)',
    structure = SMILES('[C][CH]C=C=[CH]'),
    E0 = (1089.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.89172,'amu*angstrom^2'), symmetry=1, barrier=(43.4943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (63.0773,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21402,0.0332244,-1.71911e-05,-7.36501e-09,7.12796e-12,131061,15.3217], Tmin=(100,'K'), Tmax=(920.935,'K')), NASAPolynomial(coeffs=[11.6622,0.00802887,-1.95572e-06,2.84609e-10,-1.9169e-14,128649,-33.1293], Tmin=(920.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1089.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJ3) + radical(C=C=CJ) + radical(Allyl_S)"""),
)

species(
    label = '[C]=[C]CC=[CH](19257)',
    structure = SMILES('[C]=[C]CC=[CH]'),
    E0 = (1131.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,240.145],'cm^-1')),
        HinderedRotor(inertia=(0.266673,'amu*angstrom^2'), symmetry=1, barrier=(10.7229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269836,'amu*angstrom^2'), symmetry=1, barrier=(10.7354,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30756,0.0399091,-4.49449e-05,2.98956e-08,-8.40625e-12,136181,20.1622], Tmin=(100,'K'), Tmax=(849.603,'K')), NASAPolynomial(coeffs=[6.6232,0.0195903,-9.07072e-06,1.74535e-09,-1.22727e-13,135448,0.0451475], Tmin=(849.603,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1131.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[C]#C[CH]C[CH](21433)',
    structure = SMILES('[C]#C[CH]C[CH]'),
    E0 = (1057.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,3025,407.5,1350,352.5,495.029,495.066,495.425,498.255,499.183],'cm^-1')),
        HinderedRotor(inertia=(0.0107822,'amu*angstrom^2'), symmetry=1, barrier=(1.88911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0109983,'amu*angstrom^2'), symmetry=1, barrier=(1.92456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0106848,'amu*angstrom^2'), symmetry=1, barrier=(1.88799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88021,0.0469023,-6.11633e-05,4.35103e-08,-1.17669e-11,127248,19.3532], Tmin=(100,'K'), Tmax=(1075.57,'K')), NASAPolynomial(coeffs=[8.32238,0.0143148,-3.68194e-06,4.22364e-10,-1.79315e-14,126362,-9.8755], Tmin=(1075.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1057.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCJ2_triplet) + radical(Acetyl)"""),
)

species(
    label = '[C]#C(5143)',
    structure = SMILES('[C]#C'),
    E0 = (551.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03853,0.0115449,-2.13263e-05,1.81932e-08,-5.41587e-12,66398,5.96675], Tmin=(100,'K'), Tmax=(1076.58,'K')), NASAPolynomial(coeffs=[4.00845,0.00206817,6.04935e-08,-1.17707e-10,1.2928e-14,66529.5,2.79656], Tmin=(1076.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.936,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = '[CH][CH][CH](21434)',
    structure = SMILES('[CH][CH][CH]'),
    E0 = (970.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,228.833,1326.66,1326.93,1327.2,1327.37,1327.76],'cm^-1')),
        HinderedRotor(inertia=(0.00386548,'amu*angstrom^2'), symmetry=1, barrier=(4.8339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128614,'amu*angstrom^2'), symmetry=1, barrier=(4.83416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26593,0.017361,-1.27463e-05,5.77902e-09,-1.23783e-12,116717,13.66], Tmin=(100,'K'), Tmax=(1015.96,'K')), NASAPolynomial(coeffs=[4.50709,0.0124743,-5.53133e-06,1.04452e-09,-7.27817e-14,116465,7.65249], Tmin=(1015.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(970.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJC) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[C]=[CH](19253)',
    structure = SMILES('[CH]=[C]C[C]=[CH]'),
    E0 = (1058.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,2750,2850,1437.5,1250,1305,750,350,1670,1700,300,440],'cm^-1')),
        HinderedRotor(inertia=(0.550256,'amu*angstrom^2'), symmetry=1, barrier=(12.6515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551422,'amu*angstrom^2'), symmetry=1, barrier=(12.6783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19039,0.04225,-5.14742e-05,3.63595e-08,-1.06192e-11,127386,20.3858], Tmin=(100,'K'), Tmax=(827.659,'K')), NASAPolynomial(coeffs=[7.12239,0.0184138,-8.27459e-06,1.5625e-09,-1.08408e-13,126570,-2.47531], Tmin=(827.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1058.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=C[C]=C(19277)',
    structure = SMILES('[CH]=C=C[C]=C'),
    E0 = (587.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,3010,987.5,1337.5,450,1655,1685,370],'cm^-1')),
        HinderedRotor(inertia=(1.7607,'amu*angstrom^2'), symmetry=1, barrier=(40.4819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63528,0.044168,-4.52114e-05,2.40401e-08,-4.84455e-12,70808.6,17.5611], Tmin=(100,'K'), Tmax=(1402.09,'K')), NASAPolynomial(coeffs=[11.4691,0.0097722,-1.62986e-06,9.23713e-11,5.78346e-16,68674.3,-30.9823], Tmin=(1402.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]C1C=C1(21435)',
    structure = SMILES('[CH]=[C]C1C=C1'),
    E0 = (835.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3769,0.0350407,-2.99559e-05,1.417e-08,-2.79118e-12,100582,15.5853], Tmin=(100,'K'), Tmax=(1194.62,'K')), NASAPolynomial(coeffs=[7.88925,0.0165835,-6.78048e-06,1.23677e-09,-8.46311e-14,99264.9,-11.9889], Tmin=(1194.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(835.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (867.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1032.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1142.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1046.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1285.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1089.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (931.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1325.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1324.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1397.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1082.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1072.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1005.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1505.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1069.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (1050.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1104.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (1059.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1398.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (876.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1300.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1315.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1209.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1556.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1241.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (890.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (870.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['C#C(582)', '[CH]=C=[CH](18734)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH]C=C=C=[CH](21421)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C[C]=C=[CH](21422)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.80481e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]C=[C][C]=C(19273)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(32660.1,'s^-1'), n=2.37785, Ea=(178.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH]C=[C][C]=[CH](21423)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]C1[CH]C1=[CH](21424)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.00692e+11,'s^-1'), n=0.347401, Ea=(198.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH][CH]C1[C]=C1(21425)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(221.251,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]=CC#C[CH2](19271)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH][CH](5289)', '[CH]=C=[CH](18734)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.59417e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(8)', '[CH][C]=C[C]=[CH](21426)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(8)', '[C]=[C]C=C[CH](21427)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]C=C=C[CH](21428)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]C=[C][CH2](19274)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.189e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 125 used for R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleH
Exact match found for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]C1[CH][C]=C1(21429)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_linear;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]C=[CH](19276)', '[C]=[CH](18830)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH]=[C]C=C=[CH](21430)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1968.26,'m^3/(mol*s)'), n=1.42233, Ea=(22.9562,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-H;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[CH](583)', '[CH]=C=[CH](18734)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.09425,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#C(582)', '[CH][C]=[CH](21256)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(92.9253,'m^3/(mol*s)'), n=1.51997, Ea=(27.4714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]=[C][C]=C[CH2](21431)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[CH](583)', '[CH][C]=[CH](21256)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.9156e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]=C1[CH]C=C1(21418)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[C][CH]C=C=[CH](21432)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]=[C]CC=[CH](19257)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.30972e+07,'s^-1'), n=1.70216, Ea=(184.173,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;XH_out] for rate rule [R3H_TS;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[C]#C[CH]C[CH](21433)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(645397,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_2;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[C]#C(5143)', '[CH][CH][CH](21434)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]C[C]=[CH](19253)'],
    products = ['[CH][CH]C=C=[CH](19252)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.45795e+10,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]=C=C[C]=C(19277)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH][CH]C=C=[CH](19252)'],
    products = ['[CH]=[C]C1C=C1(21435)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '4446',
    isomers = [
        '[CH][CH]C=C=[CH](19252)',
    ],
    reactants = [
        ('C#C(582)', '[CH]=C=[CH](18734)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4446',
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

