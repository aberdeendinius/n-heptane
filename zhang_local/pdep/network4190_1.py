species(
    label = '[CH2][C]=CCC=C[O](18235)',
    structure = SMILES('[CH2][C]=CCC=C[O]'),
    E0 = (374.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,310.609,310.611,310.617],'cm^-1')),
        HinderedRotor(inertia=(0.312743,'amu*angstrom^2'), symmetry=1, barrier=(21.413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312766,'amu*angstrom^2'), symmetry=1, barrier=(21.4129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312781,'amu*angstrom^2'), symmetry=1, barrier=(21.4129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916057,0.0547193,-1.47185e-05,-2.63891e-08,1.57382e-11,45193.1,27.7221], Tmin=(100,'K'), Tmax=(964.167,'K')), NASAPolynomial(coeffs=[17.308,0.017892,-5.92865e-06,1.07099e-09,-7.80429e-14,40583,-58.2765], Tmin=(964.167,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P) + radical(C=COJ)"""),
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
    label = 'C#C[CH2](17441)',
    structure = SMILES('C#C[CH2]'),
    E0 = (328.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2175,525,1131.03,1132.16,1135.9],'cm^-1')),
        HinderedRotor(inertia=(0.154206,'amu*angstrom^2'), symmetry=1, barrier=(3.5455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32026,0.0108736,8.62061e-06,-1.82973e-08,7.68649e-12,39535.3,8.27851], Tmin=(100,'K'), Tmax=(960.555,'K')), NASAPolynomial(coeffs=[6.38511,0.00814486,-2.78734e-06,4.95348e-10,-3.50148e-14,38483.6,-8.79383], Tmin=(960.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Propargyl)"""),
)

species(
    label = 'C=C1[CH]CC1[CH][O](20234)',
    structure = SMILES('C=C1[CH]CC1[CH][O]'),
    E0 = (455.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[12.1001,-0.0231216,0.000138022,-1.26988e-07,2.85982e-11,54403.5,-22.6766], Tmin=(100,'K'), Tmax=(1746.92,'K')), NASAPolynomial(coeffs=[77.4078,0.024037,-7.1366e-05,1.73744e-08,-1.28538e-12,1572.89,-460.084], Tmin=(1746.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C1C[C]=CC1(20235)',
    structure = SMILES('[O][CH]C1C[C]=CC1'),
    E0 = (492.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54571,0.0468794,-1.99043e-05,-5.79853e-10,1.72844e-12,59371.2,24.3489], Tmin=(100,'K'), Tmax=(1280.09,'K')), NASAPolynomial(coeffs=[10.7227,0.0293265,-1.237e-05,2.28436e-09,-1.56632e-13,56110.4,-25.7503], Tmin=(1280.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(492.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CCsJOH) + radical(CCOJ) + radical(cyclopentene-vinyl)"""),
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
    label = '[CH2]C#CCC=C[O](20236)',
    structure = SMILES('[CH2]C#CCC=C[O]'),
    E0 = (305.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2100,2250,500,550,3000,3100,440,815,1455,1000,396.753,397.063,398.465],'cm^-1')),
        HinderedRotor(inertia=(0.197897,'amu*angstrom^2'), symmetry=1, barrier=(22.0357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199021,'amu*angstrom^2'), symmetry=1, barrier=(22.0087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07247,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26411,0.0477099,-6.82366e-06,-3.07211e-08,1.66946e-11,36795.9,25.5563], Tmin=(100,'K'), Tmax=(959.551,'K')), NASAPolynomial(coeffs=[15.8991,0.0167667,-5.44984e-06,9.76853e-10,-7.11644e-14,32603.2,-51.6567], Tmin=(959.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=COJ) + radical(Propargyl)"""),
)

species(
    label = '[CH2][C]=CCC=C=O(20237)',
    structure = SMILES('[CH2][C]=CCC=C=O'),
    E0 = (355.581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,253.769,253.771],'cm^-1')),
        HinderedRotor(inertia=(0.204887,'amu*angstrom^2'), symmetry=1, barrier=(9.3633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.515764,'amu*angstrom^2'), symmetry=1, barrier=(23.5685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204865,'amu*angstrom^2'), symmetry=1, barrier=(9.3632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30371,0.0629454,-6.81731e-05,4.30089e-08,-1.1426e-11,42860.5,24.3366], Tmin=(100,'K'), Tmax=(897.988,'K')), NASAPolynomial(coeffs=[8.61315,0.0303856,-1.37841e-05,2.62973e-09,-1.84191e-13,41547.8,-10.1406], Tmin=(897.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = '[CH2]C=[C]CC=C[O](20238)',
    structure = SMILES('[CH2]C=[C]CC=C[O]'),
    E0 = (374.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,310.609,310.611,310.617],'cm^-1')),
        HinderedRotor(inertia=(0.312743,'amu*angstrom^2'), symmetry=1, barrier=(21.413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312766,'amu*angstrom^2'), symmetry=1, barrier=(21.4129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312781,'amu*angstrom^2'), symmetry=1, barrier=(21.4129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916057,0.0547193,-1.47185e-05,-2.63891e-08,1.57382e-11,45193.1,27.7221], Tmin=(100,'K'), Tmax=(964.167,'K')), NASAPolynomial(coeffs=[17.308,0.017892,-5.92865e-06,1.07099e-09,-7.80429e-14,40583,-58.2765], Tmin=(964.167,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CCC=[C]O(20239)',
    structure = SMILES('[CH2][C]=CCC=[C]O'),
    E0 = (473.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,326.095,326.095],'cm^-1')),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703691,0.063579,-4.5489e-05,6.8396e-09,3.96561e-12,57017.5,30.2025], Tmin=(100,'K'), Tmax=(966.336,'K')), NASAPolynomial(coeffs=[16.491,0.0182112,-6.08287e-06,1.05155e-09,-7.2807e-14,53033.3,-50.2489], Tmin=(966.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C=C[O](20240)',
    structure = SMILES('[CH2]C=C[CH]C=C[O]'),
    E0 = (238.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27925,0.0375947,4.81725e-05,-1.00985e-07,4.47848e-11,28817.6,25.6155], Tmin=(100,'K'), Tmax=(930.858,'K')), NASAPolynomial(coeffs=[20.6022,0.0114585,-1.39532e-06,1.78025e-10,-1.93609e-14,22755.2,-79.4628], Tmin=(930.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJC=C) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]=CC[C]=CO(20241)',
    structure = SMILES('[CH2][C]=CC[C]=CO'),
    E0 = (471.114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.859409,'amu*angstrom^2'), symmetry=1, barrier=(19.7595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860473,'amu*angstrom^2'), symmetry=1, barrier=(19.784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861134,'amu*angstrom^2'), symmetry=1, barrier=(19.7992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860123,'amu*angstrom^2'), symmetry=1, barrier=(19.7759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479077,0.0662053,-4.26455e-05,-3.6233e-09,9.60913e-12,56799,28.354], Tmin=(100,'K'), Tmax=(936.92,'K')), NASAPolynomial(coeffs=[19.3344,0.0137671,-3.61804e-06,5.73436e-10,-4.04467e-14,52034.2,-67.9562], Tmin=(936.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=[C]CC=C[O](20242)',
    structure = SMILES('C[C]=[C]CC=C[O]'),
    E0 = (461.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,236.68,236.68,236.68],'cm^-1')),
        HinderedRotor(inertia=(0.368295,'amu*angstrom^2'), symmetry=1, barrier=(14.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368295,'amu*angstrom^2'), symmetry=1, barrier=(14.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368295,'amu*angstrom^2'), symmetry=1, barrier=(14.6402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.691603,0.0639506,-5.49312e-05,2.42507e-08,-4.25406e-12,55581.3,28.7286], Tmin=(100,'K'), Tmax=(1375.93,'K')), NASAPolynomial(coeffs=[15.6231,0.0205428,-7.60928e-06,1.32222e-09,-8.80619e-14,51472.4,-48.0723], Tmin=(1375.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC[C]=C[O](20243)',
    structure = SMILES('[CH2]C=CC[C]=C[O]'),
    E0 = (374.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,310.609,310.611,310.617],'cm^-1')),
        HinderedRotor(inertia=(0.312743,'amu*angstrom^2'), symmetry=1, barrier=(21.413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312766,'amu*angstrom^2'), symmetry=1, barrier=(21.4129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312781,'amu*angstrom^2'), symmetry=1, barrier=(21.4129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916057,0.0547193,-1.47185e-05,-2.63891e-08,1.57382e-11,45193.1,27.7221], Tmin=(100,'K'), Tmax=(964.167,'K')), NASAPolynomial(coeffs=[17.308,0.017892,-5.92865e-06,1.07099e-09,-7.80429e-14,40583,-58.2765], Tmin=(964.167,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]=C[CH]C=CO(20244)',
    structure = SMILES('[CH2][C]=C[CH]C=CO'),
    E0 = (334.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835296,0.0491492,2.00865e-05,-7.81632e-08,3.87156e-11,40423.9,26.2732], Tmin=(100,'K'), Tmax=(915.836,'K')), NASAPolynomial(coeffs=[22.7235,0.00717471,1.00595e-06,-3.40771e-10,1.99864e-14,34165.8,-89.6781], Tmin=(915.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P) + radical(C=CCJC=C)"""),
)

species(
    label = 'C[C]=C[CH]C=C[O](20245)',
    structure = SMILES('C[C]=C[CH]C=C[O]'),
    E0 = (324.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24369,0.044739,1.45429e-05,-5.77514e-08,2.74222e-11,39197.5,25.9356], Tmin=(100,'K'), Tmax=(941.557,'K')), NASAPolynomial(coeffs=[17.1746,0.0169183,-4.63435e-06,7.87019e-10,-5.84195e-14,34430.8,-59.3448], Tmin=(941.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=CC[CH][C]=O(20246)',
    structure = SMILES('[CH2]C=CC[CH][C]=O'),
    E0 = (319.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907923,0.0595326,-4.70531e-05,1.94918e-08,-3.25348e-12,38588.9,28.1278], Tmin=(100,'K'), Tmax=(1429.39,'K')), NASAPolynomial(coeffs=[14.097,0.0226236,-8.32022e-06,1.42646e-09,-9.3788e-14,34818.5,-40.2133], Tmin=(1429.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=[C]CC=CO(20247)',
    structure = SMILES('[CH2][C]=[C]CC=CO'),
    E0 = (471.114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.859409,'amu*angstrom^2'), symmetry=1, barrier=(19.7595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860473,'amu*angstrom^2'), symmetry=1, barrier=(19.784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861134,'amu*angstrom^2'), symmetry=1, barrier=(19.7992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860123,'amu*angstrom^2'), symmetry=1, barrier=(19.7759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479077,0.0662053,-4.26455e-05,-3.6233e-09,9.60913e-12,56799,28.354], Tmin=(100,'K'), Tmax=(936.92,'K')), NASAPolynomial(coeffs=[19.3344,0.0137671,-3.61804e-06,5.73436e-10,-4.04467e-14,52034.2,-67.9562], Tmin=(936.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=CC[C]=C[O](20248)',
    structure = SMILES('C[C]=CC[C]=C[O]'),
    E0 = (461.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,236.68,236.68,236.68],'cm^-1')),
        HinderedRotor(inertia=(0.368295,'amu*angstrom^2'), symmetry=1, barrier=(14.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368295,'amu*angstrom^2'), symmetry=1, barrier=(14.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368295,'amu*angstrom^2'), symmetry=1, barrier=(14.6402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.691603,0.0639506,-5.49312e-05,2.42507e-08,-4.25406e-12,55581.3,28.7286], Tmin=(100,'K'), Tmax=(1375.93,'K')), NASAPolynomial(coeffs=[15.6231,0.0205428,-7.60928e-06,1.32222e-09,-8.80619e-14,51472.4,-48.0723], Tmin=(1375.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=CC[CH][C]=O(20249)',
    structure = SMILES('C[C]=CC[CH][C]=O'),
    E0 = (406.205,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3010,987.5,1337.5,450,1655,345.324,345.325],'cm^-1')),
        HinderedRotor(inertia=(0.00141365,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0957149,'amu*angstrom^2'), symmetry=1, barrier=(8.09965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095718,'amu*angstrom^2'), symmetry=1, barrier=(8.09963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0957198,'amu*angstrom^2'), symmetry=1, barrier=(8.09965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37421,0.0607494,-5.99697e-05,3.60017e-08,-9.30554e-12,48947.1,26.649], Tmin=(100,'K'), Tmax=(913.586,'K')), NASAPolynomial(coeffs=[7.72851,0.0329288,-1.42928e-05,2.67107e-09,-1.84954e-13,47786.1,-3.43276], Tmin=(913.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Cds_S) + radical(CCJCHO)"""),
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
    label = '[CH2][C]=C[CH2](18777)',
    structure = SMILES('[CH2][C]=C[CH2]'),
    E0 = (510.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.997629,'amu*angstrom^2'), symmetry=1, barrier=(87.0906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0766176,'amu*angstrom^2'), symmetry=1, barrier=(26.8011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.62092,0.0235622,5.21e-06,-2.26288e-08,1.0095e-11,61507.3,14.9566], Tmin=(100,'K'), Tmax=(985.296,'K')), NASAPolynomial(coeffs=[8.77433,0.0145496,-5.37933e-06,9.84612e-10,-6.99436e-14,59519.6,-18.5722], Tmin=(985.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=C(18825)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,228.264,228.889,229.07],'cm^-1')),
        HinderedRotor(inertia=(1.35219,'amu*angstrom^2'), symmetry=1, barrier=(50.6528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27541,0.0127954,9.49515e-06,-1.56026e-08,5.42938e-12,73954,11.3502], Tmin=(100,'K'), Tmax=(1063.31,'K')), NASAPolynomial(coeffs=[4.18965,0.0168435,-6.77763e-06,1.22218e-09,-8.33556e-14,73336.3,4.89309], Tmin=(1063.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C[CH]C=C[O](20250)',
    structure = SMILES('[CH2][C]=C[CH]C=C[O]'),
    E0 = (476.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,443.776,443.8,443.808],'cm^-1')),
        HinderedRotor(inertia=(0.244957,'amu*angstrom^2'), symmetry=1, barrier=(34.2094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244855,'amu*angstrom^2'), symmetry=1, barrier=(34.2102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244739,'amu*angstrom^2'), symmetry=1, barrier=(34.2088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25135,0.0421377,2.35603e-05,-7.18927e-08,3.38652e-11,57420.8,26.1893], Tmin=(100,'K'), Tmax=(930.766,'K')), NASAPolynomial(coeffs=[19.3909,0.0109953,-1.6939e-06,2.31934e-10,-2.1093e-14,52016.3,-70.9149], Tmin=(930.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CC[C]=C[O](20251)',
    structure = SMILES('[CH2][C]=CC[C]=C[O]'),
    E0 = (612.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,300.643,300.716,301.01],'cm^-1')),
        HinderedRotor(inertia=(0.282556,'amu*angstrom^2'), symmetry=1, barrier=(18.1534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282047,'amu*angstrom^2'), symmetry=1, barrier=(18.1497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85692,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.87418,0.0594249,-3.9883e-05,3.37991e-09,4.55338e-12,73796.8,28.3461], Tmin=(100,'K'), Tmax=(985.89,'K')), NASAPolynomial(coeffs=[16.1803,0.0172877,-6.14637e-06,1.10588e-09,-7.82025e-14,69808.6,-50.2], Tmin=(985.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Cds_S) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]=[C]CC=C[O](20252)',
    structure = SMILES('[CH2][C]=[C]CC=C[O]'),
    E0 = (612.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,300.643,300.716,301.01],'cm^-1')),
        HinderedRotor(inertia=(0.282556,'amu*angstrom^2'), symmetry=1, barrier=(18.1534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282047,'amu*angstrom^2'), symmetry=1, barrier=(18.1497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85692,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.87418,0.0594249,-3.9883e-05,3.37991e-09,4.55338e-12,73796.8,28.3461], Tmin=(100,'K'), Tmax=(985.89,'K')), NASAPolynomial(coeffs=[16.1803,0.0172877,-6.14637e-06,1.10588e-09,-7.82025e-14,69808.6,-50.2], Tmin=(985.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC[CH][C]=O(20253)',
    structure = SMILES('[CH2][C]=CC[CH][C]=O'),
    E0 = (557.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,285.456,4000],'cm^-1')),
        HinderedRotor(inertia=(0.351932,'amu*angstrom^2'), symmetry=1, barrier=(20.3506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5292,'amu*angstrom^2'), symmetry=1, barrier=(88.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00206943,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52981,'amu*angstrom^2'), symmetry=1, barrier=(88.4442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27781,0.0594976,-5.62055e-05,2.93436e-08,-6.31326e-12,67174.6,27.2673], Tmin=(100,'K'), Tmax=(1108.12,'K')), NASAPolynomial(coeffs=[10.5949,0.0258656,-1.06798e-05,1.95442e-09,-1.34056e-13,65109.8,-18.6386], Tmin=(1108.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC1[CH]O1(20254)',
    structure = SMILES('[CH2][C]=CCC1[CH]O1'),
    E0 = (513.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0562841,0.069892,-6.46806e-05,3.12956e-08,-5.7123e-12,61962,29.0479], Tmin=(100,'K'), Tmax=(1587.07,'K')), NASAPolynomial(coeffs=[15.7732,0.0171259,-2.64554e-06,1.27483e-10,2.28336e-15,58558.3,-49.5253], Tmin=(1587.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C1[CH]C[CH]C1[O](20255)',
    structure = SMILES('C=C1[CH]C[CH]C1[O]'),
    E0 = (391.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20697,0.0142051,0.000100932,-1.42083e-07,5.51183e-11,47198.5,22.5769], Tmin=(100,'K'), Tmax=(961.637,'K')), NASAPolynomial(coeffs=[16.4077,0.0190088,-6.19229e-06,1.25286e-09,-1.018e-13,41514,-60.7336], Tmin=(961.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclopentane) + radical(CC(C)OJ) + radical(Allyl_S) + radical(CCJCO)"""),
)

species(
    label = '[O]C1[CH]CC=[C]C1(20072)',
    structure = SMILES('[O]C1[CH]CC=[C]C1'),
    E0 = (474.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96702,0.0268269,5.45915e-05,-8.83882e-08,3.50032e-11,57142.2,22.9498], Tmin=(100,'K'), Tmax=(974.168,'K')), NASAPolynomial(coeffs=[13.6514,0.0233169,-8.47292e-06,1.62612e-09,-1.21851e-13,52755.7,-43.9445], Tmin=(974.168,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CCJCO) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=CCC=C=O(20256)',
    structure = SMILES('[CH2]C=CCC=C=O'),
    E0 = (117.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23481,0.059716,-4.88266e-05,2.14387e-08,-3.91704e-12,14261.1,24.0972], Tmin=(100,'K'), Tmax=(1275.05,'K')), NASAPolynomial(coeffs=[11.1768,0.0285271,-1.21357e-05,2.25488e-09,-1.557e-13,11725.8,-26.2829], Tmin=(1275.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=C[CH]C[CH][O](20257)',
    structure = SMILES('[CH2][C]=C[CH]C[CH][O]'),
    E0 = (700.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,3000,3050,390,425,1340,1360,335,370,324.826,324.829,2329.78,2329.78],'cm^-1')),
        HinderedRotor(inertia=(0.497739,'amu*angstrom^2'), symmetry=1, barrier=(37.2672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00967542,'amu*angstrom^2'), symmetry=1, barrier=(37.2675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.49774,'amu*angstrom^2'), symmetry=1, barrier=(37.267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144305,'amu*angstrom^2'), symmetry=1, barrier=(10.8048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92176,0.0566288,-9.91622e-06,-1.11888e-07,1.2266e-10,84290.1,24.5838], Tmin=(100,'K'), Tmax=(443.468,'K')), NASAPolynomial(coeffs=[6.11668,0.0389698,-1.84378e-05,3.52635e-09,-2.44765e-13,83719.6,5.51959], Tmin=(443.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCsJOH) + radical(Allyl_S) + radical(Allyl_P) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]=[C]CC[CH][O](20258)',
    structure = SMILES('[CH2][C]=[C]CC[CH][O]'),
    E0 = (797.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,303.534,303.558,1845.36,1845.39],'cm^-1')),
        HinderedRotor(inertia=(0.125098,'amu*angstrom^2'), symmetry=1, barrier=(8.17965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125095,'amu*angstrom^2'), symmetry=1, barrier=(8.17869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12508,'amu*angstrom^2'), symmetry=1, barrier=(8.1787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.731263,'amu*angstrom^2'), symmetry=1, barrier=(47.8043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.710386,0.0825186,-0.000131891,1.22466e-07,-4.4334e-11,95969.3,30.2759], Tmin=(100,'K'), Tmax=(843.842,'K')), NASAPolynomial(coeffs=[5.01024,0.0407458,-1.96133e-05,3.72275e-09,-2.54908e-13,96005.2,14.7742], Tmin=(843.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(797.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCsJOH) + radical(CCOJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C[CH]C[O](20259)',
    structure = SMILES('[CH2][C]=[C]C[CH]C[O]'),
    E0 = (816.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,180,994.039,994.343],'cm^-1')),
        HinderedRotor(inertia=(0.00514954,'amu*angstrom^2'), symmetry=1, barrier=(3.61337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157266,'amu*angstrom^2'), symmetry=1, barrier=(3.61585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157234,'amu*angstrom^2'), symmetry=1, barrier=(3.61512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157186,'amu*angstrom^2'), symmetry=1, barrier=(3.61401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17337,0.06857,-9.20941e-05,8.08702e-08,-2.95067e-11,98314.1,31.5624], Tmin=(100,'K'), Tmax=(796.139,'K')), NASAPolynomial(coeffs=[4.90621,0.0397402,-1.87935e-05,3.59473e-09,-2.49535e-13,98039,16.4099], Tmin=(796.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(816.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCJCO) + radical(Cds_S) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1=CCC=CO1(20171)',
    structure = SMILES('[CH2]C1=CCC=CO1'),
    E0 = (92.0067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60725,0.0341566,4.28068e-05,-8.42071e-08,3.58703e-11,11168.8,22.1507], Tmin=(100,'K'), Tmax=(949.524,'K')), NASAPolynomial(coeffs=[16.3411,0.018377,-5.38957e-06,9.72719e-10,-7.42251e-14,6284.12,-59.1564], Tmin=(949.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.0067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(1,4-Cyclohexadiene) + radical(C=C(O)CJ)"""),
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
    label = '[CH]=CCC=[C][CH2](19649)',
    structure = SMILES('[CH]=CCC=[C][CH2]'),
    E0 = (689.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0655264,'amu*angstrom^2'), symmetry=1, barrier=(20.7142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309858,'amu*angstrom^2'), symmetry=1, barrier=(7.12424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.90052,'amu*angstrom^2'), symmetry=1, barrier=(20.7047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38903,0.0502897,-3.27191e-05,8.29642e-09,-2.30657e-14,82986.9,25.1635], Tmin=(100,'K'), Tmax=(1150.67,'K')), NASAPolynomial(coeffs=[11.9984,0.022766,-9.03699e-06,1.6426e-09,-1.12829e-13,79925.9,-30.201], Tmin=(1150.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CCC=C[O](20033)',
    structure = SMILES('[CH][C]=CCC=C[O]'),
    E0 = (593.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,249.093,249.094,249.094,249.094,249.094,249.094],'cm^-1')),
        HinderedRotor(inertia=(1.11937,'amu*angstrom^2'), symmetry=1, barrier=(49.2866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11937,'amu*angstrom^2'), symmetry=1, barrier=(49.2866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11938,'amu*angstrom^2'), symmetry=1, barrier=(49.2866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.906982,0.0572433,-2.32289e-05,-1.32051e-08,9.83933e-12,71553.1,28.5831], Tmin=(100,'K'), Tmax=(987.292,'K')), NASAPolynomial(coeffs=[15.0225,0.0239288,-8.8864e-06,1.603e-09,-1.12363e-13,67602.3,-45.2282], Tmin=(987.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C1CC1[CH][O](20260)',
    structure = SMILES('C=[C]C1CC1[CH][O]'),
    E0 = (565.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21876,0.054961,-3.81222e-05,1.34978e-08,-1.9489e-12,68121.2,26.4831], Tmin=(100,'K'), Tmax=(1601.02,'K')), NASAPolynomial(coeffs=[13.2447,0.0249152,-9.97231e-06,1.7762e-09,-1.18567e-13,64270.4,-37.1952], Tmin=(1601.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCsJOH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=C=C[CH]C=C[O](20261)',
    structure = SMILES('C=C=C[CH]C=C[O]'),
    E0 = (263.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61492,'amu*angstrom^2'), symmetry=1, barrier=(37.1301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61291,'amu*angstrom^2'), symmetry=1, barrier=(37.0841,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3233,0.039683,3.0721e-05,-7.93297e-08,3.64771e-11,31832.7,24.2304], Tmin=(100,'K'), Tmax=(931.617,'K')), NASAPolynomial(coeffs=[19.4903,0.0106915,-1.51183e-06,2.05831e-10,-2.00481e-14,26320.9,-73.543], Tmin=(931.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CCJC=C)"""),
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
    label = 'C=[C]C=C(18778)',
    structure = SMILES('C=[C]C=C'),
    E0 = (292.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.58595,'amu*angstrom^2'), symmetry=1, barrier=(36.4641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60774,0.0234239,8.49718e-06,-3.01392e-08,1.43429e-11,35286.5,12.5896], Tmin=(100,'K'), Tmax=(918.282,'K')), NASAPolynomial(coeffs=[9.56656,0.0118206,-3.11021e-06,4.74837e-10,-3.20206e-14,33219.7,-24.6845], Tmin=(918.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C[CH]C=C[O](18234)',
    structure = SMILES('C=[C]C[CH]C=C[O]'),
    E0 = (375.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,478.969,478.998,479.06,479.099],'cm^-1')),
        HinderedRotor(inertia=(0.123998,'amu*angstrom^2'), symmetry=1, barrier=(20.1927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124021,'amu*angstrom^2'), symmetry=1, barrier=(20.193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123985,'amu*angstrom^2'), symmetry=1, barrier=(20.1929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05388,0.0520422,-1.07847e-05,-2.7442e-08,1.5296e-11,45258.6,26.2518], Tmin=(100,'K'), Tmax=(975.568,'K')), NASAPolynomial(coeffs=[16.2019,0.0197483,-6.97467e-06,1.28229e-09,-9.29764e-14,40884.2,-53.7259], Tmin=(975.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]C=CCC=C[O](20024)',
    structure = SMILES('[CH]C=CCC=C[O]'),
    E0 = (356.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948542,0.0525426,1.91254e-06,-4.29343e-08,2.10022e-11,42949.3,27.9603], Tmin=(100,'K'), Tmax=(967.5,'K')), NASAPolynomial(coeffs=[16.1465,0.0245391,-8.67195e-06,1.56886e-09,-1.12264e-13,38378.4,-53.2837], Tmin=(967.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]CC[C]=C[O](18236)',
    structure = SMILES('C=[C]CC[C]=C[O]'),
    E0 = (472.054,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,354.783,354.783,354.784,354.784],'cm^-1')),
        HinderedRotor(inertia=(0.173854,'amu*angstrom^2'), symmetry=1, barrier=(15.5288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173854,'amu*angstrom^2'), symmetry=1, barrier=(15.5288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173854,'amu*angstrom^2'), symmetry=1, barrier=(15.5288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.783835,0.0623133,-4.60197e-05,1.21238e-08,6.52725e-13,56898,28.8019], Tmin=(100,'K'), Tmax=(1035.07,'K')), NASAPolynomial(coeffs=[15.2131,0.0213279,-8.03774e-06,1.45223e-09,-1.00928e-13,53119.4,-45.132], Tmin=(1035.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(472.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CCC=C[O](18238)',
    structure = SMILES('[CH]=[C]CCC=C[O]'),
    E0 = (481.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,278.981,278.987,278.987],'cm^-1')),
        HinderedRotor(inertia=(0.309373,'amu*angstrom^2'), symmetry=1, barrier=(17.0869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309339,'amu*angstrom^2'), symmetry=1, barrier=(17.0869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309364,'amu*angstrom^2'), symmetry=1, barrier=(17.0868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3854.44,'J/mol'), sigma=(6.42147,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.05 K, Pc=33.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769055,0.0609423,-3.69303e-05,-8.86212e-10,6.12015e-12,58013.3,28.7837], Tmin=(100,'K'), Tmax=(987.657,'K')), NASAPolynomial(coeffs=[16.4338,0.0193423,-6.92307e-06,1.25043e-09,-8.85111e-14,53853.7,-51.9882], Tmin=(987.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(481.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]CC[CH][C]=O(18237)',
    structure = SMILES('C=[C]CC[CH][C]=O'),
    E0 = (418.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1855,455,950,3025,407.5,1350,352.5,342.448,342.843,342.957],'cm^-1')),
        HinderedRotor(inertia=(0.15818,'amu*angstrom^2'), symmetry=1, barrier=(13.1979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00143696,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158214,'amu*angstrom^2'), symmetry=1, barrier=(13.1979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00956198,'amu*angstrom^2'), symmetry=1, barrier=(79.1642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30612,0.060838,-5.66798e-05,3.053e-08,-6.95933e-12,50423.6,27.3274], Tmin=(100,'K'), Tmax=(1034.59,'K')), NASAPolynomial(coeffs=[9.08844,0.0307499,-1.30573e-05,2.42108e-09,-1.67159e-13,48813.3,-10.4826], Tmin=(1034.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CCC=CO(20031)',
    structure = SMILES('[CH][C]=CCC=CO'),
    E0 = (452.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509865,0.0640524,-2.61181e-05,-2e-08,1.47838e-11,54555.4,28.598], Tmin=(100,'K'), Tmax=(943.829,'K')), NASAPolynomial(coeffs=[18.1655,0.0204258,-6.36769e-06,1.07274e-09,-7.47827e-14,49832.9,-62.9214], Tmin=(943.829,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C1C[CH]C1[O](20262)',
    structure = SMILES('C=[C]C1C[CH]C1[O]'),
    E0 = (581.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76826,0.0311469,4.56947e-05,-8.28832e-08,3.41047e-11,70065.8,26.3428], Tmin=(100,'K'), Tmax=(964.977,'K')), NASAPolynomial(coeffs=[15.0821,0.0204908,-6.96385e-06,1.32015e-09,-9.98969e-14,65422.9,-48.1575], Tmin=(964.977,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C=C[CH]C=CO(20263)',
    structure = SMILES('C=C=C[CH]C=CO'),
    E0 = (122.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.906533,0.0467041,2.72076e-05,-8.55384e-08,4.12957e-11,14835.8,24.3168], Tmin=(100,'K'), Tmax=(917.215,'K')), NASAPolynomial(coeffs=[22.8223,0.00687185,1.18754e-06,-3.66776e-10,2.10239e-14,8470.73,-92.3023], Tmin=(917.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=[C]CCC=C=O(18226)',
    structure = SMILES('C=[C]CCC=C=O'),
    E0 = (216.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2120,512.5,787.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,288.937,289.108,289.21],'cm^-1')),
        HinderedRotor(inertia=(0.00202029,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196954,'amu*angstrom^2'), symmetry=1, barrier=(11.7078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196714,'amu*angstrom^2'), symmetry=1, barrier=(11.7087,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54353,0.0561705,-4.40564e-05,1.90855e-08,-3.56638e-12,26111,25.0685], Tmin=(100,'K'), Tmax=(1212.47,'K')), NASAPolynomial(coeffs=[8.84295,0.0320894,-1.42646e-05,2.70474e-09,-1.88828e-13,24341,-11.5533], Tmin=(1212.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([C]=C)C=C[O](20264)',
    structure = SMILES('[CH2]C([C]=C)C=C[O]'),
    E0 = (430.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,323.318,323.405,323.77],'cm^-1')),
        HinderedRotor(inertia=(0.18517,'amu*angstrom^2'), symmetry=1, barrier=(13.7618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347635,'amu*angstrom^2'), symmetry=1, barrier=(25.7913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185193,'amu*angstrom^2'), symmetry=1, barrier=(13.7313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.663469,0.0675075,-6.04857e-05,2.52744e-08,-3.05725e-12,51935.5,27.9098], Tmin=(100,'K'), Tmax=(970.67,'K')), NASAPolynomial(coeffs=[14.7366,0.0212903,-7.26322e-06,1.21925e-09,-8.08434e-14,48648.7,-42.4237], Tmin=(970.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C]C1CC=CO1(20181)',
    structure = SMILES('C=[C]C1CC=CO1'),
    E0 = (176.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72644,0.0248127,8.14462e-05,-1.3757e-07,5.92264e-11,21340.7,19.8083], Tmin=(100,'K'), Tmax=(912.135,'K')), NASAPolynomial(coeffs=[20.6795,0.00739601,2.04715e-06,-5.73005e-10,3.4936e-14,15150.1,-84.8676], Tmin=(912.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[CH]CC=O(20130)',
    structure = SMILES('[CH2][C]=C[CH]CC=O'),
    E0 = (430.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,361.045,3257.82],'cm^-1')),
        HinderedRotor(inertia=(0.258144,'amu*angstrom^2'), symmetry=1, barrier=(23.8787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.087239,'amu*angstrom^2'), symmetry=1, barrier=(8.06979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0872385,'amu*angstrom^2'), symmetry=1, barrier=(8.06979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0399789,'amu*angstrom^2'), symmetry=1, barrier=(38.6483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58281,0.0553308,-4.17944e-05,1.72187e-08,-3.06481e-12,51816,28.0862], Tmin=(100,'K'), Tmax=(1260.11,'K')), NASAPolynomial(coeffs=[8.8305,0.0323242,-1.4408e-05,2.72977e-09,-1.90271e-13,49989.4,-8.55537], Tmin=(1260.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CCC[C]=O(18801)',
    structure = SMILES('[CH2][C]=CCC[C]=O'),
    E0 = (390.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1855,455,950,3000,3100,440,815,1455,1000,292.677,292.831],'cm^-1')),
        HinderedRotor(inertia=(0.115808,'amu*angstrom^2'), symmetry=1, barrier=(7.03406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115439,'amu*angstrom^2'), symmetry=1, barrier=(7.03265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00196504,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266467,'amu*angstrom^2'), symmetry=1, barrier=(16.2199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27284,0.0638107,-6.57225e-05,4.04753e-08,-1.06761e-11,47022.1,27.5359], Tmin=(100,'K'), Tmax=(897.357,'K')), NASAPolynomial(coeffs=[8.06496,0.0335346,-1.51138e-05,2.87698e-09,-2.01322e-13,45803.1,-4.49653], Tmin=(897.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=[C]CCC=O(20265)',
    structure = SMILES('[CH2][C]=[C]CCC=O'),
    E0 = (468.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,308.964,2292.93],'cm^-1')),
        HinderedRotor(inertia=(0.140506,'amu*angstrom^2'), symmetry=1, barrier=(9.51882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140462,'amu*angstrom^2'), symmetry=1, barrier=(9.51879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140548,'amu*angstrom^2'), symmetry=1, barrier=(9.5198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.626573,'amu*angstrom^2'), symmetry=1, barrier=(42.4363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74791,0.0581748,-3.55741e-05,-3.07411e-08,4.60604e-11,56366.5,25.5991], Tmin=(100,'K'), Tmax=(497.396,'K')), NASAPolynomial(coeffs=[6.04447,0.0379582,-1.78391e-05,3.43359e-09,-2.40613e-13,55761.8,6.08868], Tmin=(497.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC[CH][CH]O1(20152)',
    structure = SMILES('[CH2]C1=CC[CH][CH]O1'),
    E0 = (356.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45305,0.0368178,4.10228e-05,-9.13129e-08,4.17362e-11,42994.9,21.9999], Tmin=(100,'K'), Tmax=(910.76,'K')), NASAPolynomial(coeffs=[18.685,0.0115824,-4.98819e-07,-1.03311e-10,5.66245e-15,37763.8,-71.0096], Tmin=(910.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(3,4-Dihydro-2H-pyran) + radical(CCJCO) + radical(C=C(O)CJ) + radical(CCsJOC(O))"""),
)

species(
    label = '[C]1=CC[CH][CH]OC1(20017)',
    structure = SMILES('[C]1=CC[CH][CH]OC1'),
    E0 = (499.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74947,-0.00134554,0.000204767,-3.0567e-07,1.29984e-10,60224.3,26.3166], Tmin=(100,'K'), Tmax=(899.534,'K')), NASAPolynomial(coeffs=[37.1463,-0.0249929,2.11615e-05,-4.29606e-09,2.85447e-13,48444.7,-170.784], Tmin=(899.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(CCJCO) + radical(CCsJOCs) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC=CC=O(20266)',
    structure = SMILES('[CH2]C=CC=CC=O'),
    E0 = (52.4633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19817,0.051949,-2.18837e-05,-5.89878e-09,4.90491e-12,6418.88,24.4996], Tmin=(100,'K'), Tmax=(1111.19,'K')), NASAPolynomial(coeffs=[13.2992,0.0262162,-1.12126e-05,2.13734e-09,-1.51482e-13,2628.95,-40.109], Tmin=(1111.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.4633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]=CC([CH2])C=O(18134)',
    structure = SMILES('[CH2][C]=CC([CH2])C=O'),
    E0 = (433.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,255.82],'cm^-1')),
        HinderedRotor(inertia=(0.132739,'amu*angstrom^2'), symmetry=1, barrier=(6.16371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132733,'amu*angstrom^2'), symmetry=1, barrier=(6.16336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132743,'amu*angstrom^2'), symmetry=1, barrier=(6.16404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.787026,'amu*angstrom^2'), symmetry=1, barrier=(36.5475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3682.09,'J/mol'), sigma=(6.19766,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=575.13 K, Pc=35.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.854116,0.0744442,-9.81694e-05,7.86217e-08,-2.62642e-11,52285.9,27.7555], Tmin=(100,'K'), Tmax=(767.354,'K')), NASAPolynomial(coeffs=[7.87421,0.0352295,-1.63905e-05,3.12222e-09,-2.16736e-13,51285.7,-3.75057], Tmin=(767.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C=C1[CH]CC1C=O(20106)',
    structure = SMILES('C=C1[CH]CC1C=O'),
    E0 = (127.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.818,-0.026257,0.000145876,-1.33215e-07,3.01879e-11,14984.8,-23.2777], Tmin=(100,'K'), Tmax=(1724.02,'K')), NASAPolynomial(coeffs=[73.226,0.0283053,-7.30323e-05,1.77279e-08,-1.31354e-12,-35471.2,-437.907], Tmin=(1724.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S)"""),
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
    label = '[CH]CC=[C][CH2](18854)',
    structure = SMILES('[CH]CC=[C][CH2]'),
    E0 = (785.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,374.118,374.461,374.553,1854.97],'cm^-1')),
        HinderedRotor(inertia=(0.00120176,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00120334,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372728,'amu*angstrom^2'), symmetry=1, barrier=(37.0812,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82458,0.0420401,-2.86433e-05,8.66561e-09,-6.57495e-13,94497.9,21.4552], Tmin=(100,'K'), Tmax=(1180.82,'K')), NASAPolynomial(coeffs=[10.592,0.0188325,-7.40916e-06,1.3331e-09,-9.07983e-14,91974.8,-24.216], Tmin=(1180.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCJ2_triplet) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CCCC=O(20047)',
    structure = SMILES('[CH][C]=CCCC=O'),
    E0 = (449.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86095,0.0550025,-1.53968e-05,-5.1432e-08,5.24604e-11,54119.3,25.5519], Tmin=(100,'K'), Tmax=(489.268,'K')), NASAPolynomial(coeffs=[4.43006,0.0454274,-2.10795e-05,4.053e-09,-2.85182e-13,53731.1,13.596], Tmin=(489.268,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C1C[CH][CH]O1(20190)',
    structure = SMILES('C=[C]C1C[CH][CH]O1'),
    E0 = (486.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73228,0.0313342,5.32437e-05,-1.02826e-07,4.63368e-11,58657.3,23.6024], Tmin=(100,'K'), Tmax=(891.847,'K')), NASAPolynomial(coeffs=[16.904,0.0131523,-4.29237e-08,-3.01997e-10,2.40005e-14,53968.1,-58.9734], Tmin=(891.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(CCsJOCs) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=CCC=O(20267)',
    structure = SMILES('C=[C]C=CCC=O'),
    E0 = (157.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14085,0.052814,-2.28167e-05,-5.61669e-09,4.8368e-12,19069.7,24.168], Tmin=(100,'K'), Tmax=(1123.45,'K')), NASAPolynomial(coeffs=[13.9442,0.0257035,-1.12871e-05,2.17957e-09,-1.55476e-13,15027,-44.2801], Tmin=(1123.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]CC=CC=O(18246)',
    structure = SMILES('C=[C]CC=CC=O'),
    E0 = (203.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,258.737,262.176],'cm^-1')),
        HinderedRotor(inertia=(0.259532,'amu*angstrom^2'), symmetry=1, barrier=(12.647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265547,'amu*angstrom^2'), symmetry=1, barrier=(12.6782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258735,'amu*angstrom^2'), symmetry=1, barrier=(12.6349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8411,0.0524533,-3.41345e-05,1.05484e-08,-1.34473e-12,24500,24.1019], Tmin=(100,'K'), Tmax=(1683.28,'K')), NASAPolynomial(coeffs=[11.0729,0.0305155,-1.45851e-05,2.80575e-09,-1.9478e-13,21392.1,-25.2433], Tmin=(1683.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1CC1C=O(20138)',
    structure = SMILES('C=[C]C1CC1C=O'),
    E0 = (239.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58891,0.0414677,5.27076e-06,-3.46815e-08,1.56279e-11,28941.2,24.2047], Tmin=(100,'K'), Tmax=(1001.41,'K')), NASAPolynomial(coeffs=[12.5246,0.0245273,-9.40912e-06,1.75674e-09,-1.2581e-13,25410.1,-35.2638], Tmin=(1001.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S)"""),
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
    E0 = (374.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (459.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (492.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (531.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (570.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (472.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (578.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (635.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (522.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (621.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (622.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (475.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (489.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (478.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (419.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (534.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (603.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (481.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (733.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (704.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (688.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (824.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (824.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (769.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (588.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (412.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (474.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (399.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (723.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (860.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (841.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (382.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (1096.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (805.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (565.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (483.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (535.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (550.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (570.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (569.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (614.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (626.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (462.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (603.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (581.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (399.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (399.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (588.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (381.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (537.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (546.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (547.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (592.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (433.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (499.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (452.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (593.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (382.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (852.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (493.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (486.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (397.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (402.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (382.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=CC=O(5269)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=C1[CH]CC1[CH][O](20234)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(85.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[O][CH]C1C[C]=CC1(20235)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.13873e+09,'s^-1'), n=0.337103, Ea=(118.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 116.5 to 118.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2]C#CCC=C[O](20236)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH2][C]=CCC=C=O(20237)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C=C[O](5266)', 'C#C[CH2](17441)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0565962,'m^3/(mol*s)'), n=2.418, Ea=(54.0165,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CsJ-CdHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C=[C]CC=C[O](20238)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]=CCC=[C]O(20239)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH2]C=C[CH]C=C[O](20240)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=CC[C]=CO(20241)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[C]=[C]CC=C[O](20242)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C=CC[C]=C[O](20243)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSD;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH2][C]=C[CH]C=CO(20244)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.26335e+06,'s^-1'), n=1.7925, Ea=(114.6,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;Y_rad_out;Cs_H_out_H/OneDe] + [R4H_SDS;O_rad_out;Cs_H_out_1H] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C[C]=C[CH]C=C[O](20245)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.13003e+06,'s^-1'), n=1.69583, Ea=(104.251,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R4HJ_1;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH2]C=CC[CH][C]=O(20246)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSD;Cd_rad_out;XH_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=[C]CC=CO(20247)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(313415,'s^-1'), n=1.7968, Ea=(63.8264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_RSMS;Cd_rad_out;XH_out] + [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[C]=CC[C]=C[O](20248)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[C]=CC[CH][C]=O(20249)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.31883e+06,'s^-1'), n=1.02765, Ea=(75.0925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6HJ_4;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C[O](602)', '[CH2][C]=C[CH2](18777)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.169e+08,'m^3/(mol*s)'), n=-0.455312, Ea=(0.377199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_pri_rad;Cd_rad] + [C_rad/H2/Cd;Y_rad] for rate rule [C_rad/H2/Cd;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C=C[O](5266)', '[CH][C]=C(18825)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cd]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2][C]=C[CH]C=C[O](20250)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.14e-13,'cm^3/(molecule*s)'), n=0.611, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 53 used for C_rad/H/CdCd;H_rad
Exact match found for rate rule [C_rad/H/CdCd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.8 to -1.8 kJ/mol.
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH2][C]=CC[C]=C[O](20251)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2][C]=[C]CC=C[O](20252)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH2][C]=CC[CH][C]=O(20253)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH2][C]=CCC1[CH]O1(20254)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=C1[CH]C[CH]C1[O](20255)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.80657e+08,'s^-1'), n=0.835, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra_pri;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[O]C1[CH]CC=[C]C1(20072)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.15136e+07,'s^-1'), n=0.706682, Ea=(99.6302,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 95.8 to 99.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH2]C=CCC=C=O(20256)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=C[CH]C[CH][O](20257)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]=[C]CC[CH][O](20258)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=[C]C[CH]C[O](20259)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH2]C1=CCC=CO1(20171)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(T)(63)', '[CH]=CCC=[C][CH2](19649)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(8)', '[CH][C]=CCC=C[O](20033)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=[C]C1CC1[CH][O](20260)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.18499e+11,'s^-1'), n=-0.0609598, Ea=(190.783,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 190.3 to 190.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(8)', 'C=C=C[CH]C=C[O](20261)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(8)', 'C#C[CH]CC=C[O](20022)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C[O](602)', 'C=[C]C=C(18778)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.123469,'m^3/(mol*s)'), n=2.00579, Ea=(36.0234,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]C[CH]C=C[O](18234)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.01225e+10,'s^-1'), n=0.845153, Ea=(195.403,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH]C=CCC=C[O](20024)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C]CC[C]=C[O](18236)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]CCC=C[O](18238)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=[C]CC[CH][C]=O(18237)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH][C]=CCC=CO(20031)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=[C]C1C[CH]C1[O](20262)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.66065e+07,'s^-1'), n=1.25778, Ea=(207.022,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 203.8 to 207.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=C=C[CH]C=CO(20263)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=[C]CCC=C=O(18226)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C([C]=C)C=C[O](20264)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=[C]C1CC=CO1(20181)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=CC=O(5269)', '[CH][C]=C(18825)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(0.0102751,'m^3/(mol*s)'), n=2.40501, Ea=(4.48561,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDeH;CJ] for rate rule [Cds-HH_Cds-COH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][C]=C[CH]CC=O(20130)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.169e+11,'s^-1'), n=0.707, Ea=(116.068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2][C]=CCC[C]=O(18801)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2][C]=[C]CCC=O(20265)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH2]C1=CC[CH][CH]O1(20152)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(58.6284,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[C]1=CC[CH][CH]OC1(20017)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2.97949e+09,'s^-1'), n=0.649948, Ea=(124.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cs2H] + [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 118.3 to 124.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['[CH2]C=CC=CC=O(20266)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]=CC([CH2])C=O(18134)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=C1[CH]CC1C=O(20106)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH]=O(373)', '[CH]CC=[C][CH2](18854)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH][C]=CCCC=O(20047)'],
    products = ['[CH2][C]=CCC=C[O](18235)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R5Hall;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=[C]C1C[CH][CH]O1(20190)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.14364e+09,'s^-1'), n=0.582609, Ea=(112.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 109.2 to 112.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=[C]C=CCC=O(20267)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=[C]CC=CC=O(18246)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2][C]=CCC=C[O](18235)'],
    products = ['C=[C]C1CC1C=O(20138)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '4190',
    isomers = [
        '[CH2][C]=CCC=C[O](18235)',
    ],
    reactants = [
        ('C=CC=O(5269)', 'C#C[CH2](17441)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4190',
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

