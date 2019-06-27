species(
    label = '[CH2]C(CO[O])C([O])=O(6933)',
    structure = SMILES('[CH2]C(CO[O])C([O])=O'),
    E0 = (1.34773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.631229,0.0871474,-0.000151736,1.47403e-07,-5.46628e-11,270.798,29.6376], Tmin=(100,'K'), Tmax=(845.837,'K')), NASAPolynomial(coeffs=[3.80397,0.0429894,-2.17252e-05,4.18196e-09,-2.87661e-13,776.971,21.0271], Tmin=(845.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.34773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CJC(C)C=O) + radical(ROOJ)"""),
)

species(
    label = 'O=C=O(1731)',
    structure = SMILES('O=C=O'),
    E0 = (-403.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.923,1087.69,1087.69,2296.71],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2028.74,'J/mol'), sigma=(3.763,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(2.65,'angstroms^3'), rotrelaxcollnum=2.1, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27862,0.0027414,7.16119e-06,-1.08033e-08,4.14308e-12,-48470.3,5.97933], Tmin=(100,'K'), Tmax=(988.876,'K')), NASAPolynomial(coeffs=[4.54605,0.0029192,-1.15488e-06,2.27663e-10,-1.70918e-14,-48980.3,-1.43251], Tmin=(988.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.087,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cdd-OdOd)"""),
)

species(
    label = 'C=CCO[O](6082)',
    structure = SMILES('C=CCO[O]'),
    E0 = (76.4976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.316477,'amu*angstrom^2'), symmetry=1, barrier=(7.27644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.315948,'amu*angstrom^2'), symmetry=1, barrier=(7.26427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3494.07,'J/mol'), sigma=(5.81539,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=545.77 K, Pc=40.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56676,0.0335761,-2.36362e-05,9.70939e-09,-1.82265e-12,9250.36,17.5875], Tmin=(100,'K'), Tmax=(1155.37,'K')), NASAPolynomial(coeffs=[5.58064,0.0231419,-1.00897e-05,1.8929e-09,-1.31328e-13,8553.92,2.61201], Tmin=(1155.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.4976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = '[O]OCC1CC1([O])[O](7509)',
    structure = SMILES('[O]OCC1CC1([O])[O]'),
    E0 = (152.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59865,0.0580246,-5.21452e-05,2.79103e-08,-6.70623e-12,18396.7,25.7707], Tmin=(100,'K'), Tmax=(948.293,'K')), NASAPolynomial(coeffs=[6.82636,0.0359734,-1.72645e-05,3.3884e-09,-2.41432e-13,17405.2,0.827646], Tmin=(948.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = '[CH2]C1COOC1([O])[O](7510)',
    structure = SMILES('[CH2]C1COOC1([O])[O]'),
    E0 = (118.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85632,0.049424,-2.4812e-05,-1.46913e-08,1.85171e-11,14344.6,25.7937], Tmin=(100,'K'), Tmax=(647.071,'K')), NASAPolynomial(coeffs=[7.27663,0.0293336,-9.34043e-06,1.41129e-09,-8.41976e-14,13362.3,-0.166998], Tmin=(647.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsOsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CCOJ) + radical(CCOJ) + radical(Isobutyl)"""),
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
    label = 'C=C(CO[O])C([O])=O(7511)',
    structure = SMILES('C=C(CO[O])C([O])=O'),
    E0 = (-1.10896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,1079.71,1600,2381.13,3200],'cm^-1')),
        HinderedRotor(inertia=(0.134195,'amu*angstrom^2'), symmetry=1, barrier=(3.0854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134195,'amu*angstrom^2'), symmetry=1, barrier=(3.0854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134195,'amu*angstrom^2'), symmetry=1, barrier=(3.0854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20993,0.0771648,-0.000147932,1.54876e-07,-5.99484e-11,-48.1656,28.7793], Tmin=(100,'K'), Tmax=(851.108,'K')), NASAPolynomial(coeffs=[0.0383994,0.0445256,-2.31806e-05,4.50041e-09,-3.10197e-13,1532.84,42.3588], Tmin=(851.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.10896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]O[O](61)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (200.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1055.29],'cm^-1')),
        HinderedRotor(inertia=(0.00752578,'amu*angstrom^2'), symmetry=1, barrier=(5.89392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3200.22,'J/mol'), sigma=(5.39124,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=499.87 K, Pc=46.34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38991,0.0147959,-1.79906e-05,1.54065e-08,-5.43781e-12,24103.1,9.19559], Tmin=(100,'K'), Tmax=(865.677,'K')), NASAPolynomial(coeffs=[3.82258,0.00980114,-4.14553e-06,7.46963e-10,-4.99013e-14,24140.4,7.8189], Tmin=(865.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsHHH) + radical(ROOJ) + radical(CsJOOH)"""),
)

species(
    label = 'C=CC([O])=O(3003)',
    structure = SMILES('C=CC([O])=O'),
    E0 = (-28.1395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180,1961.5,1961.57,1961.7],'cm^-1')),
        HinderedRotor(inertia=(0.148974,'amu*angstrom^2'), symmetry=1, barrier=(3.4252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89654,0.0323031,-5.7836e-05,6.80187e-08,-2.92995e-11,-3352.53,16.0958], Tmin=(100,'K'), Tmax=(822.704,'K')), NASAPolynomial(coeffs=[-0.923254,0.0300379,-1.57145e-05,3.10029e-09,-2.17304e-13,-2018.85,38.0643], Tmin=(822.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.1395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[O][C]=O(2059)',
    structure = SMILES('[O][C]=O'),
    E0 = (33.3014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81048,-0.00025715,1.76446e-05,-2.38747e-08,9.15883e-12,4016.03,8.55818], Tmin=(100,'K'), Tmax=(975.962,'K')), NASAPolynomial(coeffs=[6.50409,-1.44217e-05,-6.90664e-08,7.0435e-11,-9.1126e-15,2952.93,-7.12421], Tmin=(975.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-OdOsH) + radical(OJC=O) + radical((O)CJOH)"""),
)

species(
    label = '[CH2][CH]CO[O](6081)',
    structure = SMILES('[CH2][CH]CO[O]'),
    E0 = (348.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.00295793,'amu*angstrom^2'), symmetry=1, barrier=(5.15327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223054,'amu*angstrom^2'), symmetry=1, barrier=(5.12845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222928,'amu*angstrom^2'), symmetry=1, barrier=(5.12556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03778,0.0476625,-6.74371e-05,5.88574e-08,-2.07062e-11,42019.6,22.7181], Tmin=(100,'K'), Tmax=(834.77,'K')), NASAPolynomial(coeffs=[5.09101,0.024729,-1.13075e-05,2.11543e-09,-1.4439e-13,41799.2,10.2722], Tmin=(834.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = 'C[C](CO[O])C([O])=O(7512)',
    structure = SMILES('C[C](CO[O])C([O])=O'),
    E0 = (-56.5891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20358,0.0703392,-0.000104369,1.00145e-07,-3.86311e-11,-6714,28.1198], Tmin=(100,'K'), Tmax=(805.695,'K')), NASAPolynomial(coeffs=[2.89502,0.0439919,-2.18986e-05,4.25411e-09,-2.96889e-13,-6403.96,23.9405], Tmin=(805.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.5891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCJ(C)CO) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([CH]OO)C([O])=O(7513)',
    structure = SMILES('[CH2]C([CH]OO)C([O])=O'),
    E0 = (37.9255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.357039,0.0936551,-0.000164568,1.59197e-07,-5.90209e-11,4679.48,30.0115], Tmin=(100,'K'), Tmax=(837.041,'K')), NASAPolynomial(coeffs=[4.65668,0.0438112,-2.27453e-05,4.42699e-09,-3.06562e-13,4986.01,16.1637], Tmin=(837.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.9255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCsJOOH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CC([CH]O[O])C([O])=O(7514)',
    structure = SMILES('CC([CH]O[O])C([O])=O'),
    E0 = (-20.5805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.899652,0.0801376,-0.000133862,1.31249e-07,-4.96247e-11,-2375.2,28.4719], Tmin=(100,'K'), Tmax=(836.309,'K')), NASAPolynomial(coeffs=[2.75312,0.0446369,-2.24143e-05,4.32458e-09,-2.9876e-13,-1753.74,25.4302], Tmin=(836.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.5805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2][C](CO[O])C(=O)O(7515)',
    structure = SMILES('[CH2][C](CO[O])C(=O)O'),
    E0 = (-71.7835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582196,0.0827309,-0.000126189,1.09981e-07,-3.82792e-11,-8517.77,30.5364], Tmin=(100,'K'), Tmax=(817.509,'K')), NASAPolynomial(coeffs=[7.89659,0.0350985,-1.70599e-05,3.26597e-09,-2.25562e-13,-9317.92,-0.856868], Tmin=(817.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.7835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCJ(C)CO) + radical(ROOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([CH]O[O])C(=O)O(7516)',
    structure = SMILES('[CH2]C([CH]O[O])C(=O)O'),
    E0 = (-35.7749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.274619,0.092586,-0.000155963,1.41618e-07,-4.96054e-11,-4178.81,30.9009], Tmin=(100,'K'), Tmax=(852.094,'K')), NASAPolynomial(coeffs=[7.73134,0.0357842,-1.75995e-05,3.34215e-09,-2.27911e-13,-4658.26,0.763434], Tmin=(852.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.7749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCsJOOH) + radical(ROOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(COO)=C([O])[O](7517)',
    structure = SMILES('[CH2]C(COO)=C([O])[O]'),
    E0 = (-46.0671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469767,0.0807196,-9.60843e-05,5.95245e-08,-1.47929e-11,-5415.99,27.8699], Tmin=(100,'K'), Tmax=(976.249,'K')), NASAPolynomial(coeffs=[13.7501,0.0263063,-1.24793e-05,2.43222e-09,-1.72662e-13,-8008.99,-35.8809], Tmin=(976.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.0671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C([CH2])C([O])=O(3081)',
    structure = SMILES('[CH2]C([CH2])C([O])=O'),
    E0 = (148.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,732.939,732.946,732.962],'cm^-1')),
        HinderedRotor(inertia=(0.00589214,'amu*angstrom^2'), symmetry=1, barrier=(2.24621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00590098,'amu*angstrom^2'), symmetry=1, barrier=(2.24954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0978921,'amu*angstrom^2'), symmetry=1, barrier=(2.25073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4058.35,'J/mol'), sigma=(6.30806,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=633.90 K, Pc=36.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51763,0.0636488,-0.000104602,1.01684e-07,-3.83728e-11,17902.1,21.8906], Tmin=(100,'K'), Tmax=(831.284,'K')), NASAPolynomial(coeffs=[3.33648,0.0351773,-1.76443e-05,3.4105e-09,-2.36121e-13,18281,17.5499], Tmin=(831.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJC(C)C=O) + radical(CJC(C)C=O) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C=C([O])[O](3060)',
    structure = SMILES('[CH2]C=C([O])[O]'),
    E0 = (74.6188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,402.923,402.926],'cm^-1')),
        HinderedRotor(inertia=(0.301493,'amu*angstrom^2'), symmetry=1, barrier=(34.7355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18883,0.0335239,-1.94151e-05,-2.53461e-09,4.26037e-12,9045.37,16.414], Tmin=(100,'K'), Tmax=(992.311,'K')), NASAPolynomial(coeffs=[11.6135,0.00923936,-3.42475e-06,6.42069e-10,-4.67977e-14,6500.12,-32.3821], Tmin=(992.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.6188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(CO[O])=C([O])[O](7518)',
    structure = SMILES('[CH2]C(CO[O])=C([O])[O]'),
    E0 = (105.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,322.226,322.449,2062.16],'cm^-1')),
        HinderedRotor(inertia=(0.188422,'amu*angstrom^2'), symmetry=1, barrier=(13.9047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188586,'amu*angstrom^2'), symmetry=1, barrier=(13.9045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.825512,'amu*angstrom^2'), symmetry=1, barrier=(60.8756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.617215,0.0785271,-0.000105823,7.56183e-08,-2.15881e-11,12859.5,28.0977], Tmin=(100,'K'), Tmax=(856.247,'K')), NASAPolynomial(coeffs=[12.1411,0.0246952,-1.15233e-05,2.20103e-09,-1.53407e-13,10885.9,-25.7107], Tmin=(856.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ) + radical(ROOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH]O[O])C([O])=O(7519)',
    structure = SMILES('[CH2]C([CH]O[O])C([O])=O'),
    E0 = (189.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612545,0.090147,-0.000169444,1.68428e-07,-6.25435e-11,22950.3,29.854], Tmin=(100,'K'), Tmax=(858.506,'K')), NASAPolynomial(coeffs=[3.08037,0.0421303,-2.17425e-05,4.18348e-09,-2.86204e-13,23872.3,26.1625], Tmin=(858.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJC(C)C=O) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]C(CO[O])[C]1OO1(7520)',
    structure = SMILES('[CH2]C(CO[O])[C]1OO1'),
    E0 = (354.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.637428,0.0691151,-6.45046e-05,2.35566e-08,-3.58912e-14,42800.5,31.5025], Tmin=(100,'K'), Tmax=(872.216,'K')), NASAPolynomial(coeffs=[15.8559,0.0167717,-4.49413e-06,6.24215e-10,-3.68832e-14,39482,-43.6418], Tmin=(872.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(dioxirane) + radical(ROOJ) + radical(Cs_P) + radical(Isobutyl)"""),
)

species(
    label = '[O]OCC1CO[C]1[O](7521)',
    structure = SMILES('[O]OCC1CO[C]1[O]'),
    E0 = (165.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52837,0.0576169,-4.70901e-05,1.04888e-08,7.32619e-12,19962.5,27.2813], Tmin=(100,'K'), Tmax=(656.525,'K')), NASAPolynomial(coeffs=[8.14013,0.0295097,-1.06916e-05,1.77742e-09,-1.13777e-13,18832,-3.83293], Tmin=(656.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Oxetane) + radical(ROOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1COOO[C]1[O](7522)',
    structure = SMILES('[CH2]C1COOO[C]1[O]'),
    E0 = (307.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.421974,0.0669273,-6.44714e-05,3.33806e-08,-6.641e-12,37105,27.5111], Tmin=(100,'K'), Tmax=(1410.62,'K')), NASAPolynomial(coeffs=[14.3167,0.0187784,-3.96875e-06,3.90063e-10,-1.4984e-14,34055.4,-41.2177], Tmin=(1410.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(123trioxane) + radical(CCOJ) + radical(Isobutyl) + radical(Cs_P)"""),
)

species(
    label = 'C=C(CO[O])C(=O)O(6930)',
    structure = SMILES('C=C(CO[O])C(=O)O'),
    E0 = (-226.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874838,0.0795691,-0.000134325,1.27902e-07,-4.69436e-11,-27177.4,29.8161], Tmin=(100,'K'), Tmax=(841.391,'K')), NASAPolynomial(coeffs=[4.67649,0.0382022,-1.9051e-05,3.66231e-09,-2.52177e-13,-26992.6,17.0317], Tmin=(841.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C=C(COO)C([O])=O(7523)',
    structure = SMILES('C=C(COO)C([O])=O'),
    E0 = (-153.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953445,0.0806908,-0.000143157,1.45857e-07,-5.65658e-11,-18318.9,28.9399], Tmin=(100,'K'), Tmax=(829.725,'K')), NASAPolynomial(coeffs=[1.59929,0.0462333,-2.4199e-05,4.74765e-09,-3.30867e-13,-17347.2,32.4464], Tmin=(829.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C](CO[O])C([O])[O](7524)',
    structure = SMILES('[CH2][C](CO[O])C([O])[O]'),
    E0 = (397.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1933.02,1933.06,1933.13,1934.32],'cm^-1')),
        HinderedRotor(inertia=(2.57289,'amu*angstrom^2'), symmetry=1, barrier=(59.1558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192936,'amu*angstrom^2'), symmetry=1, barrier=(5.1161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22293,'amu*angstrom^2'), symmetry=1, barrier=(5.1256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222276,'amu*angstrom^2'), symmetry=1, barrier=(5.11056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21946,0.0753696,-0.000134101,1.37414e-07,-5.26145e-11,47916.5,37.0435], Tmin=(100,'K'), Tmax=(858.3,'K')), NASAPolynomial(coeffs=[0.201339,0.0458463,-2.26166e-05,4.30415e-09,-2.93754e-13,49353.5,49.1528], Tmin=(858.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(C2CJCOOH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([CH]O[O])C([O])[O](7525)',
    structure = SMILES('[CH2]C([CH]O[O])C([O])[O]'),
    E0 = (395.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,180,1412.39,1412.72,1412.83],'cm^-1')),
        HinderedRotor(inertia=(0.162579,'amu*angstrom^2'), symmetry=1, barrier=(3.73801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162264,'amu*angstrom^2'), symmetry=1, barrier=(3.73077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162299,'amu*angstrom^2'), symmetry=1, barrier=(3.73156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162403,'amu*angstrom^2'), symmetry=1, barrier=(3.73396,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.629497,0.0914743,-0.000173823,1.74694e-07,-6.48096e-11,47654.8,34.4664], Tmin=(100,'K'), Tmax=(873.138,'K')), NASAPolynomial(coeffs=[1.58411,0.0452444,-2.24956e-05,4.248e-09,-2.8675e-13,49083.6,39.1271], Tmin=(873.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCOJ) + radical(ROOJ) + radical(Isobutyl) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]C([CH]O[O])[C]([O])O(7526)',
    structure = SMILES('[CH2]C([CH]O[O])[C]([O])O'),
    E0 = (374.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,222.863,822.863,1782.9],'cm^-1')),
        HinderedRotor(inertia=(0.122458,'amu*angstrom^2'), symmetry=1, barrier=(3.49607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122458,'amu*angstrom^2'), symmetry=1, barrier=(3.49607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122458,'amu*angstrom^2'), symmetry=1, barrier=(3.49607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122458,'amu*angstrom^2'), symmetry=1, barrier=(3.49607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122458,'amu*angstrom^2'), symmetry=1, barrier=(3.49607,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314508,0.0948225,-0.000171436,1.60995e-07,-5.66779e-11,45209,36.1183], Tmin=(100,'K'), Tmax=(881.579,'K')), NASAPolynomial(coeffs=[6.08288,0.0365086,-1.75276e-05,3.25042e-09,-2.167e-13,45440.9,16.1], Tmin=(881.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_P) + radical(CCOJ) + radical(CCsJOOH) + radical(Isobutyl)"""),
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
    label = '[O]C(=O)C1COC1(7527)',
    structure = SMILES('[O]C(=O)C1COC1'),
    E0 = (-248.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20101,0.035955,-6.16132e-06,-1.14869e-08,6.14262e-12,-29836.6,18.8507], Tmin=(100,'K'), Tmax=(933.124,'K')), NASAPolynomial(coeffs=[5.6552,0.0299356,-1.06111e-05,1.78452e-09,-1.16929e-13,-30863.8,0.375348], Tmin=(933.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Oxetane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1COOC1=O(7528)',
    structure = SMILES('[CH2]C1COOC1=O'),
    E0 = (-203.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95644,0.0236679,6.55296e-05,-1.0625e-07,4.31418e-11,-24425,22.0797], Tmin=(100,'K'), Tmax=(957.288,'K')), NASAPolynomial(coeffs=[16.8233,0.0149896,-4.61272e-06,9.16167e-10,-7.49318e-14,-29720.1,-61.7855], Tmin=(957.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclopentane) + radical(CJC(C)C=O)"""),
)

species(
    label = '[O]OCCC=C([O])[O](7529)',
    structure = SMILES('[O]OCCC=C([O])[O]'),
    E0 = (-35.8437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676144,0.0775796,-9.90501e-05,6.91062e-08,-1.95502e-11,-4195.26,28.7122], Tmin=(100,'K'), Tmax=(858.866,'K')), NASAPolynomial(coeffs=[11.2159,0.0284923,-1.33188e-05,2.55932e-09,-1.79393e-13,-6005.69,-20.5323], Tmin=(858.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.8437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]OC[CH]CC([O])=O(6932)',
    structure = SMILES('[O]OC[CH]CC([O])=O'),
    E0 = (-2.98751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931478,0.0812662,-0.000142385,1.43777e-07,-5.51199e-11,-262.162,30.6503], Tmin=(100,'K'), Tmax=(840.914,'K')), NASAPolynomial(coeffs=[1.47637,0.0467351,-2.38177e-05,4.61122e-09,-3.18611e-13,775.464,34.8305], Tmin=(840.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.98751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = '[O]O(16)',
    structure = SMILES('[O]O'),
    E0 = (-8.19602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1036.72,2034.11,2034.11],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.04595,-0.00173474,1.0377e-05,-1.02207e-08,3.3493e-12,-986.755,4.63579], Tmin=(100,'K'), Tmax=(932.129,'K')), NASAPolynomial(coeffs=[3.21022,0.00367946,-1.27704e-06,2.18051e-10,-1.46343e-14,-910.359,8.18305], Tmin=(932.129,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.19602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsH) + group(O2s-OsH) + radical(HOOJ)"""),
)

species(
    label = '[CH2]C(=C)C([O])=O(6553)',
    structure = SMILES('[CH2]C(=C)C([O])=O'),
    E0 = (84.8686,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,1343.27,1343.71,1344.45,1344.75],'cm^-1')),
        HinderedRotor(inertia=(0.135472,'amu*angstrom^2'), symmetry=1, barrier=(3.11477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135221,'amu*angstrom^2'), symmetry=1, barrier=(3.109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.87128,0.0287993,-1.40179e-05,2.55385e-09,-1.57424e-13,10183.9,12.5909], Tmin=(100,'K'), Tmax=(3036.44,'K')), NASAPolynomial(coeffs=[40.7824,-0.0129792,3.23888e-06,-4.92499e-10,3.2261e-14,-15387.6,-211.677], Tmin=(3036.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.8686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C=O)CJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OCC1COC1=O(7505)',
    structure = SMILES('[O]OCC1COC1=O'),
    E0 = (-262.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55131,0.0455824,-8.62083e-06,-2.31528e-08,1.34019e-11,-31454.4,25.932], Tmin=(100,'K'), Tmax=(920.498,'K')), NASAPolynomial(coeffs=[11.609,0.0233166,-7.27501e-06,1.17577e-09,-7.77109e-14,-34214.4,-26.6914], Tmin=(920.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-262.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(ROOJ)"""),
)

species(
    label = '[O]C(=O)C1COOC1(7530)',
    structure = SMILES('[O]C(=O)C1COOC1'),
    E0 = (-276.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61883,0.0456454,-1.10723e-05,-1.62852e-08,9.72849e-12,-33192.5,22.4195], Tmin=(100,'K'), Tmax=(936.22,'K')), NASAPolynomial(coeffs=[9.6025,0.0281334,-9.60853e-06,1.60936e-09,-1.06659e-13,-35414.8,-19.4558], Tmin=(936.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(12dioxolane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1COOOC1=O(7531)',
    structure = SMILES('[CH2]C1COOOC1=O'),
    E0 = (-189.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37995,0.0296316,7.69765e-05,-1.27291e-07,5.10985e-11,-22648.9,23.8524], Tmin=(100,'K'), Tmax=(981.215,'K')), NASAPolynomial(coeffs=[22.2504,0.0142056,-5.92317e-06,1.38011e-09,-1.18073e-13,-30097.7,-93.5264], Tmin=(981.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C=C([O])OCO[O](7532)',
    structure = SMILES('[CH2]C=C([O])OCO[O]'),
    E0 = (-16.6067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03272,'amu*angstrom^2'), symmetry=1, barrier=(23.7443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03179,'amu*angstrom^2'), symmetry=1, barrier=(23.723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03186,'amu*angstrom^2'), symmetry=1, barrier=(23.7244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03247,'amu*angstrom^2'), symmetry=1, barrier=(23.7386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66931,0.0912774,-9.89188e-05,4.82282e-08,-8.51901e-12,-1764.34,33.9482], Tmin=(100,'K'), Tmax=(1660.89,'K')), NASAPolynomial(coeffs=[27.3085,-0.00036303,3.57944e-06,-8.35094e-10,5.84358e-14,-8376.18,-111.482], Tmin=(1660.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.6067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C[O])C([O])=O(7533)',
    structure = SMILES('[CH2]C(C[O])C([O])=O'),
    E0 = (3.54318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40511,0.0710858,-0.000128357,1.35237e-07,-5.33075e-11,505.998,26.1984], Tmin=(100,'K'), Tmax=(841.667,'K')), NASAPolynomial(coeffs=[-0.486256,0.0462335,-2.3755e-05,4.61335e-09,-3.19208e-13,2023.03,42.1179], Tmin=(841.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.54318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJC(C)C=O) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O]OCC=C([O])[O](7534)',
    structure = SMILES('[O]OCC=C([O])[O]'),
    E0 = (-6.50653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,324.729,324.75,324.759],'cm^-1')),
        HinderedRotor(inertia=(0.103376,'amu*angstrom^2'), symmetry=1, barrier=(7.73628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103377,'amu*angstrom^2'), symmetry=1, barrier=(7.73658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34191,0.0625221,-8.9579e-05,7.01863e-08,-2.21385e-11,-690.634,24.4532], Tmin=(100,'K'), Tmax=(775.211,'K')), NASAPolynomial(coeffs=[9.28741,0.0215237,-1.02476e-05,1.96182e-09,-1.3625e-13,-1922.5,-11.856], Tmin=(775.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.50653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([C]=O)CO[O](7535)',
    structure = SMILES('[CH2]C([C]=O)CO[O]'),
    E0 = (200.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1390,370,380,2900,435,1855,455,950,2750,2850,1437.5,1250,1305,750,350,213.472],'cm^-1')),
        HinderedRotor(inertia=(0.21563,'amu*angstrom^2'), symmetry=1, barrier=(7.06615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21931,'amu*angstrom^2'), symmetry=1, barrier=(7.05532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220434,'amu*angstrom^2'), symmetry=1, barrier=(7.06806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217939,'amu*angstrom^2'), symmetry=1, barrier=(7.05162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796839,0.0785509,-0.000128662,1.13789e-07,-3.89448e-11,24164.1,28.7269], Tmin=(100,'K'), Tmax=(859.465,'K')), NASAPolynomial(coeffs=[7.94107,0.0289991,-1.37292e-05,2.56842e-09,-1.7369e-13,23538.2,-1.15516], Tmin=(859.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(ROOJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH]C(CO[O])C([O])=O(7536)',
    structure = SMILES('[CH]C(CO[O])C([O])=O'),
    E0 = (239.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848392,0.0803882,-0.00013626,1.3045e-07,-4.82133e-11,28854.1,29.5111], Tmin=(100,'K'), Tmax=(836.886,'K')), NASAPolynomial(coeffs=[4.49813,0.039201,-1.98814e-05,3.84233e-09,-2.65398e-13,29074.7,17.5206], Tmin=(836.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(ROOJ) + radical(CCJ2_triplet)"""),
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
    E0 = (1.34773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (152.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (118.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (221.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (203.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (123.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1.34773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (146.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (167.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (118.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (112.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (110.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (131.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (151.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (274.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (382.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (317.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (401.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (354.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (165.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (307.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (64.7479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (26.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (419.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (458.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (383.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (84.2458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (94.5464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (161.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (124.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (140.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (9.63205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (8.46053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (10.2178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (297.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (251.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (409.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (606.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (450.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['O=C=O(1731)', 'C=CCO[O](6082)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[O]OCC1CC1([O])[O](7509)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.93521e+09,'s^-1'), n=0.743095, Ea=(150.928,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 147.9 to 150.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[CH2]C1COOC1([O])[O](7510)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.14284e+07,'s^-1'), n=0.933356, Ea=(117.294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.73205080757
family: Intra_R_Add_Exocyclic
Ea raised from 113.1 to 117.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C=C(CO[O])C([O])=O(7511)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(72.1434,'m^3/(mol*s)'), n=1.66666, Ea=(10.8177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;HJ] for rate rule [Cds-COCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]O[O](61)', 'C=CC([O])=O(3003)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00178011,'m^3/(mol*s)'), n=2.4143, Ea=(31.0889,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeH_Cds;CsJ] for rate rule [Cds-COH_Cds;CsJ-OsHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][C]=O(2059)', 'C=CCO[O](6082)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.540096,'m^3/(mol*s)'), n=2.05449, Ea=(13.6169,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-HH;CJ] + [Cds-Cs\O2s/H_Cds-HH;YJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C=O(1731)', '[CH2][CH]CO[O](6081)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(55.615,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CJ] for rate rule [CO2;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 54.2 to 55.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[C](CO[O])C([O])=O(7512)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[CH2]C([CH]OO)C([O])=O(7513)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.66e+08,'s^-1'), n=1.28, Ea=(166.272,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 241 used for R3H_SS_O;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['CC([CH]O[O])C([O])=O(7514)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 339 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[CH2][C](CO[O])C(=O)O(7515)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.99823e+07,'s^-1'), n=1.57622, Ea=(111.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[CH2]C([CH]O[O])C(=O)O(7516)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.44e+08,'s^-1'), n=1.09, Ea=(109.37,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[CH2]C(COO)=C([O])[O](7517)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.17136e+07,'s^-1'), n=1.54267, Ea=(130.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS_OCs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O2(2)', '[CH2]C([CH2])C([O])=O(3081)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.66237e+08,'m^3/(mol*s)'), n=-0.783071, Ea=(11.596,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for O2_birad;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]O[O](61)', '[CH2]C=C([O])[O](3060)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.88087e+06,'m^3/(mol*s)'), n=0.114385, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_pri_rad] for rate rule [Y_rad;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C]=O(2059)', '[CH2][CH]CO[O](6081)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH2]C(CO[O])=C([O])[O](7518)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH2]C([CH]O[O])C([O])=O(7519)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[CH2]C(CO[O])[C]1OO1(7520)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.29014e+11,'s^-1'), n=0.514092, Ea=(353.47,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 351.6 to 353.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[O]OCC1CO[C]1[O](7521)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.54148e+08,'s^-1'), n=0.924088, Ea=(163.914,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 160.7 to 163.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[CH2]C1COOO[C]1[O](7522)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(927160,'s^-1'), n=1.14062, Ea=(306.006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 303.0 to 306.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['C=C(CO[O])C(=O)O(6930)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['C=C(COO)C([O])=O(7523)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][C](CO[O])C([O])[O](7524)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH]O[O])C([O])[O](7525)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH]O[O])[C]([O])O(7526)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['O(T)(63)', '[O]C(=O)C1COC1(7527)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_pri_rad_intra;OO] for rate rule [R3OO_SS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['O(T)(63)', '[CH2]C1COOC1=O(7528)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.82826e+10,'s^-1'), n=0.54, Ea=(93.1987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R4OO;Y_rad_intra;OO] for rate rule [R4OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[O]OCCC=C([O])[O](7529)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[O]OC[CH]CC([O])=O(6932)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[O]O(16)', '[CH2]C(=C)C([O])=O(6553)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[O]OCC1COC1=O(7505)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[O]C(=O)C1COOC1(7530)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(CO[O])C([O])=O(6933)'],
    products = ['[CH2]C1COOOC1=O(7531)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.42e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C=C([O])OCO[O](7532)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(T)(63)', '[CH2]C(C[O])C([O])=O(7533)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CH2(T)(28)', '[O]OCC=C([O])[O](7534)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O(T)(63)', '[CH2]C([C]=O)CO[O](7535)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(8)', '[CH]C(CO[O])C([O])=O(7536)'],
    products = ['[CH2]C(CO[O])C([O])=O(6933)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '2032',
    isomers = [
        '[CH2]C(CO[O])C([O])=O(6933)',
    ],
    reactants = [
        ('O=C=O(1731)', 'C=CCO[O](6082)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '2032',
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

