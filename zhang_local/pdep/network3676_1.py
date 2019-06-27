species(
    label = '[CH2]C(C=O)CC=C[O](12914)',
    structure = SMILES('[CH2]C(C=O)CC=C[O]'),
    E0 = (-9.46149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,560.661,560.713],'cm^-1')),
        HinderedRotor(inertia=(0.212843,'amu*angstrom^2'), symmetry=1, barrier=(14.8561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646117,'amu*angstrom^2'), symmetry=1, barrier=(14.8555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646125,'amu*angstrom^2'), symmetry=1, barrier=(14.8557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0665881,'amu*angstrom^2'), symmetry=1, barrier=(14.8556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00791887,0.0803703,-7.48883e-05,3.57969e-08,-6.79269e-12,-988.013,32.4133], Tmin=(100,'K'), Tmax=(1275.49,'K')), NASAPolynomial(coeffs=[17.6727,0.0249735,-9.74143e-06,1.74671e-09,-1.18836e-13,-5494.32,-57.1077], Tmin=(1275.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.46149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJC(C)C=O)"""),
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
    label = '[O]C=CCC1CC1[O](15400)',
    structure = SMILES('[O]C=CCC1CC1[O]'),
    E0 = (80.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.693858,0.0498516,3.29604e-05,-8.95457e-08,4.08576e-11,9765.69,28.4495], Tmin=(100,'K'), Tmax=(944.522,'K')), NASAPolynomial(coeffs=[22.3803,0.0154438,-3.6062e-06,6.41971e-10,-5.34805e-14,3107.15,-88.4991], Tmin=(944.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O][CH]C1CC(C=O)C1(15401)',
    structure = SMILES('[O][CH]C1CC(C=O)C1'),
    E0 = (114.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10869,0.0597657,-3.56821e-05,1.01995e-08,-1.17314e-12,13924.4,28.0864], Tmin=(100,'K'), Tmax=(1944.1,'K')), NASAPolynomial(coeffs=[15.7932,0.0295525,-1.23708e-05,2.20571e-09,-1.45197e-13,8214.7,-52.5204], Tmin=(1944.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1CC=COC1[O](15402)',
    structure = SMILES('[CH2]C1CC=COC1[O]'),
    E0 = (43.3482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3281,0.0398779,4.46803e-05,-9.16329e-08,4.00783e-11,5327.22,24.1638], Tmin=(100,'K'), Tmax=(926.947,'K')), NASAPolynomial(coeffs=[16.5496,0.0225857,-5.64606e-06,8.82397e-10,-6.30274e-14,426.317,-59.3305], Tmin=(926.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.3482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(CCOJ) + radical(Isobutyl)"""),
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
    label = 'C=C(C=O)CC=C[O](15403)',
    structure = SMILES('C=C(C=O)CC=C[O]'),
    E0 = (-109.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,245.384,245.413,245.45],'cm^-1')),
        HinderedRotor(inertia=(0.443866,'amu*angstrom^2'), symmetry=1, barrier=(18.9692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.443848,'amu*angstrom^2'), symmetry=1, barrier=(18.969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.443849,'amu*angstrom^2'), symmetry=1, barrier=(18.9692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.461114,0.0660566,-4.03183e-05,3.00153e-09,3.62076e-12,-13018.7,29.4262], Tmin=(100,'K'), Tmax=(1089.77,'K')), NASAPolynomial(coeffs=[17.7149,0.0234293,-1.01405e-05,1.97285e-09,-1.42394e-13,-18008.5,-60.937], Tmin=(1089.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C=O)CC=C=O(15404)',
    structure = SMILES('[CH2]C(C=O)CC=C=O'),
    E0 = (-27.2931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5,207.397,207.399],'cm^-1')),
        HinderedRotor(inertia=(0.211071,'amu*angstrom^2'), symmetry=1, barrier=(6.44217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211036,'amu*angstrom^2'), symmetry=1, barrier=(6.44167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00391977,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499164,'amu*angstrom^2'), symmetry=1, barrier=(15.2374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.838707,0.0768969,-8.70372e-05,5.07166e-08,-8.13215e-12,-3175.62,28.941], Tmin=(100,'K'), Tmax=(612.073,'K')), NASAPolynomial(coeffs=[8.35092,0.0382679,-1.80145e-05,3.46991e-09,-2.43309e-13,-4291.24,-5.21459], Tmin=(612.073,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.2931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + radical(CJC(C)C=O)"""),
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
    label = 'C=CCC=C[O](12633)',
    structure = SMILES('C=CCC=C[O]'),
    E0 = (21.4197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.92568,'amu*angstrom^2'), symmetry=1, barrier=(21.2832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.927956,'amu*angstrom^2'), symmetry=1, barrier=(21.3355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61321,0.0368322,2.28043e-05,-6.16688e-08,2.78452e-11,2676.41,22.6597], Tmin=(100,'K'), Tmax=(948.598,'K')), NASAPolynomial(coeffs=[16.1169,0.0147806,-4.16463e-06,7.44334e-10,-5.71701e-14,-1834.7,-55.8207], Tmin=(948.598,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.4197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'CC(=C[O])CC=C[O](15405)',
    structure = SMILES('CC(=C[O])CC=C[O]'),
    E0 = (-84.9653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.959466,'amu*angstrom^2'), symmetry=1, barrier=(22.06,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958357,'amu*angstrom^2'), symmetry=1, barrier=(22.0345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960483,'amu*angstrom^2'), symmetry=1, barrier=(22.0834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.152985,0.0655224,-8.5205e-06,-4.98996e-08,2.77895e-11,-10063,29.927], Tmin=(100,'K'), Tmax=(939.923,'K')), NASAPolynomial(coeffs=[23.4536,0.0142626,-3.15886e-06,5.16512e-10,-4.13221e-14,-16559,-92.297], Tmin=(939.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.9653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C=O)CC=[C]O(15406)',
    structure = SMILES('[CH2]C(C=O)CC=[C]O'),
    E0 = (88.8201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,352.114,352.172],'cm^-1')),
        HinderedRotor(inertia=(0.136571,'amu*angstrom^2'), symmetry=1, barrier=(12.0251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136579,'amu*angstrom^2'), symmetry=1, barrier=(12.0246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136661,'amu*angstrom^2'), symmetry=1, barrier=(12.0248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136537,'amu*angstrom^2'), symmetry=1, barrier=(12.0243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136553,'amu*angstrom^2'), symmetry=1, barrier=(12.0239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0548331,0.0861689,-9.5014e-05,5.54452e-08,-1.29155e-11,10825.2,33.9648], Tmin=(100,'K'), Tmax=(1044.53,'K')), NASAPolynomial(coeffs=[15.2593,0.0279436,-1.1399e-05,2.07809e-09,-1.42435e-13,7648.95,-40.0504], Tmin=(1044.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.8201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CJC(C)C=O) + radical(C=CJO)"""),
)

species(
    label = 'CC([CH]C=C[O])C=O(15407)',
    structure = SMILES('CC([CH]C=C[O])C=O'),
    E0 = (-20.0702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.485171,0.0652447,-3.12649e-05,-9.73424e-09,9.32843e-12,-2276.68,32.8426], Tmin=(100,'K'), Tmax=(996.564,'K')), NASAPolynomial(coeffs=[17.1672,0.0243745,-9.01529e-06,1.64983e-09,-1.1711e-13,-6897.07,-54.0808], Tmin=(996.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CCJCC=O) + radical(C=COJ)"""),
)

species(
    label = 'CC([C]=O)CC=C[O](15408)',
    structure = SMILES('CC([C]=O)CC=C[O]'),
    E0 = (-61.2836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.348652,0.0691176,-4.17416e-05,6.42669e-10,5.77497e-12,-7229.41,31.9929], Tmin=(100,'K'), Tmax=(1005.66,'K')), NASAPolynomial(coeffs=[17.4086,0.0242168,-9.00849e-06,1.64098e-09,-1.15661e-13,-11821.5,-56.1792], Tmin=(1005.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.2836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CC(C)CJ=O) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C=O)C[C]=CO(15409)',
    structure = SMILES('[CH2]C(C=O)C[C]=CO'),
    E0 = (86.9177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,235.71,235.711],'cm^-1')),
        HinderedRotor(inertia=(0.374215,'amu*angstrom^2'), symmetry=1, barrier=(14.7538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374211,'amu*angstrom^2'), symmetry=1, barrier=(14.7538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374212,'amu*angstrom^2'), symmetry=1, barrier=(14.7538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374213,'amu*angstrom^2'), symmetry=1, barrier=(14.7538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374213,'amu*angstrom^2'), symmetry=1, barrier=(14.7538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.379679,0.0913321,-0.000101309,5.7176e-08,-1.26125e-11,10615.8,32.8648], Tmin=(100,'K'), Tmax=(1113.55,'K')), NASAPolynomial(coeffs=[18.866,0.0221988,-8.18304e-06,1.42239e-09,-9.53349e-14,6329.61,-62.0545], Tmin=(1113.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.9177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CJC(C)C=O) + radical(Cds_S)"""),
)

species(
    label = 'CC(C=O)C[C]=C[O](15410)',
    structure = SMILES('CC(C=O)C[C]=C[O]'),
    E0 = (17.8696,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,3010,987.5,1337.5,450,1655,268.922,270.641,272.119],'cm^-1')),
        HinderedRotor(inertia=(0.0023406,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305284,'amu*angstrom^2'), symmetry=1, barrier=(15.6433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305123,'amu*angstrom^2'), symmetry=1, barrier=(15.6178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.302913,'amu*angstrom^2'), symmetry=1, barrier=(15.6339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215659,0.0754702,-6.566e-05,2.9514e-08,-5.30808e-12,2291.89,31.7917], Tmin=(100,'K'), Tmax=(1335.23,'K')), NASAPolynomial(coeffs=[16.6952,0.0261019,-1.01996e-05,1.82327e-09,-1.23455e-13,-2108.93,-52.477], Tmin=(1335.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.8696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=O)C=C[CH]O(15411)',
    structure = SMILES('[CH2]C(C=O)C=C[CH]O'),
    E0 = (8.64778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.640279,0.0747072,-6.84509e-05,3.38948e-08,-6.95421e-12,1160.36,32.1209], Tmin=(100,'K'), Tmax=(1149.93,'K')), NASAPolynomial(coeffs=[12.4559,0.0336067,-1.48379e-05,2.81282e-09,-1.96816e-13,-1557.06,-26.5334], Tmin=(1149.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.64778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(C=CCJO)"""),
)

species(
    label = 'CC(C=O)C[CH][C]=O(15412)',
    structure = SMILES('CC(C=O)C[CH][C]=O'),
    E0 = (-35.7338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85744,0.0725716,-7.13665e-05,4.17062e-08,-1.04214e-11,-4187.55,29.8905], Tmin=(100,'K'), Tmax=(945.12,'K')), NASAPolynomial(coeffs=[8.99692,0.038121,-1.66867e-05,3.13384e-09,-2.17732e-13,-5726.01,-8.9179], Tmin=(945.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.7338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = '[CH2]C(=C[O])CC=CO(15413)',
    structure = SMILES('[CH2]C(=C[O])CC=CO'),
    E0 = (-74.9287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254689,0.0699241,-2.94804e-06,-7.03467e-08,3.90964e-11,-8836.64,30.2619], Tmin=(100,'K'), Tmax=(917.918,'K')), NASAPolynomial(coeffs=[28.9985,0.00452591,2.47746e-06,-6.10335e-10,3.70053e-14,-16822.3,-122.608], Tmin=(917.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.9287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)CC=CO(15414)',
    structure = SMILES('[CH2]C([C]=O)CC=CO'),
    E0 = (7.76453,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,234.704,235.298],'cm^-1')),
        HinderedRotor(inertia=(0.396904,'amu*angstrom^2'), symmetry=1, barrier=(15.6001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.397838,'amu*angstrom^2'), symmetry=1, barrier=(15.5975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398226,'amu*angstrom^2'), symmetry=1, barrier=(15.5992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.397879,'amu*angstrom^2'), symmetry=1, barrier=(15.5978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398997,'amu*angstrom^2'), symmetry=1, barrier=(15.5987,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.674437,0.0898926,-9.38849e-05,4.8588e-08,-9.66492e-12,1113.15,34.6083], Tmin=(100,'K'), Tmax=(1293.49,'K')), NASAPolynomial(coeffs=[22.1765,0.0160702,-4.61433e-06,6.903e-10,-4.26688e-14,-4534.16,-80.4937], Tmin=(1293.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.76453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2][CH]CC=C[O](6404)',
    structure = SMILES('[CH2][CH]CC=C[O]'),
    E0 = (292.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180,484.096,1122.21],'cm^-1')),
        HinderedRotor(inertia=(0.0140785,'amu*angstrom^2'), symmetry=1, barrier=(2.43178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0151486,'amu*angstrom^2'), symmetry=1, barrier=(2.58464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94874,'amu*angstrom^2'), symmetry=1, barrier=(21.8134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5584,0.0454766,-1.65477e-05,-1.12403e-08,7.76691e-12,35222,26.8772], Tmin=(100,'K'), Tmax=(989.998,'K')), NASAPolynomial(coeffs=[11.9567,0.021597,-7.8421e-06,1.3995e-09,-9.72263e-14,32274.4,-27.6723], Tmin=(989.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ) + radical(C=COJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(=C[O])CC=C[O](15415)',
    structure = SMILES('[CH2]C(=C[O])CC=C[O]'),
    E0 = (66.5339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20374,'amu*angstrom^2'), symmetry=1, barrier=(27.6763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20284,'amu*angstrom^2'), symmetry=1, barrier=(27.6556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20297,'amu*angstrom^2'), symmetry=1, barrier=(27.6587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.160755,0.0629193,5.04292e-07,-6.40534e-08,3.42392e-11,8160.27,30.1802], Tmin=(100,'K'), Tmax=(930.93,'K')), NASAPolynomial(coeffs=[25.6705,0.00833852,-2.17803e-07,-3.87095e-11,-3.98446e-15,1026.21,-103.871], Tmin=(930.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.5339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(C=COJ)"""),
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
    label = '[CH2]C([CH2])C=O(6124)',
    structure = SMILES('[CH2]C([CH2])C=O'),
    E0 = (188.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.235764,'amu*angstrom^2'), symmetry=1, barrier=(5.42068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235141,'amu*angstrom^2'), symmetry=1, barrier=(5.40636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.526268,'amu*angstrom^2'), symmetry=1, barrier=(12.0999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64923,0.0560988,-7.86058e-05,6.57621e-08,-2.23429e-11,22710.5,19.9459], Tmin=(100,'K'), Tmax=(819.334,'K')), NASAPolynomial(coeffs=[6.60175,0.0257593,-1.17819e-05,2.21158e-09,-1.51531e-13,22105.8,-1.6983], Tmin=(819.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([CH]C=C[O])C=O(15416)',
    structure = SMILES('[CH2]C([CH]C=C[O])C=O'),
    E0 = (190.441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,579.073,579.134],'cm^-1')),
        HinderedRotor(inertia=(0.00911797,'amu*angstrom^2'), symmetry=1, barrier=(2.16985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0675106,'amu*angstrom^2'), symmetry=1, barrier=(16.0672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357991,'amu*angstrom^2'), symmetry=1, barrier=(16.0672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.698814,'amu*angstrom^2'), symmetry=1, barrier=(16.0671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0329734,0.077113,-7.29327e-05,3.47641e-08,-6.49332e-12,23056.1,34.8227], Tmin=(100,'K'), Tmax=(1307.04,'K')), NASAPolynomial(coeffs=[18.7387,0.0198671,-7.23584e-06,1.25493e-09,-8.39865e-14,18166.2,-60.4305], Tmin=(1307.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCJCC=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=O)C[C]=C[O](15417)',
    structure = SMILES('[CH2]C(C=O)C[C]=C[O]'),
    E0 = (228.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,386.239,386.554,387.761],'cm^-1')),
        HinderedRotor(inertia=(0.128629,'amu*angstrom^2'), symmetry=1, barrier=(13.6413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128672,'amu*angstrom^2'), symmetry=1, barrier=(13.6418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128807,'amu*angstrom^2'), symmetry=1, barrier=(13.6409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128516,'amu*angstrom^2'), symmetry=1, barrier=(13.6401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280626,0.081363,-8.7125e-05,4.90069e-08,-1.10435e-11,27602.2,31.9098], Tmin=(100,'K'), Tmax=(1075.16,'K')), NASAPolynomial(coeffs=[14.8378,0.0272053,-1.15678e-05,2.15703e-09,-1.49853e-13,24471.9,-39.3749], Tmin=(1075.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CJC(C)C=O) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([C]=O)CC=C[O](15418)',
    structure = SMILES('[CH2]C([C]=O)CC=C[O]'),
    E0 = (149.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,390.356,390.356,390.356],'cm^-1')),
        HinderedRotor(inertia=(0.136724,'amu*angstrom^2'), symmetry=1, barrier=(14.7842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136725,'amu*angstrom^2'), symmetry=1, barrier=(14.7842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136726,'amu*angstrom^2'), symmetry=1, barrier=(14.7842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136725,'amu*angstrom^2'), symmetry=1, barrier=(14.7842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0242105,0.0795362,-7.86039e-05,3.92808e-08,-7.70329e-12,18097.7,33.5109], Tmin=(100,'K'), Tmax=(1245.05,'K')), NASAPolynomial(coeffs=[18.3095,0.0207904,-7.82849e-06,1.38362e-09,-9.36933e-14,13544.5,-58.7125], Tmin=(1245.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C(C=O)C[CH][C]=O(15419)',
    structure = SMILES('[CH2]C(C=O)C[CH][C]=O'),
    E0 = (174.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,311.101,2494.98],'cm^-1')),
        HinderedRotor(inertia=(0.130754,'amu*angstrom^2'), symmetry=1, barrier=(8.98955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00174067,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130762,'amu*angstrom^2'), symmetry=1, barrier=(8.98929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130894,'amu*angstrom^2'), symmetry=1, barrier=(8.98986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.917312,'amu*angstrom^2'), symmetry=1, barrier=(63.0672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.818485,0.0784602,-8.49431e-05,3.37737e-08,7.47242e-12,21127.3,30.4469], Tmin=(100,'K'), Tmax=(579.327,'K')), NASAPolynomial(coeffs=[9.3989,0.0355198,-1.59754e-05,2.98622e-09,-2.04821e-13,19859.5,-8.62596], Tmin=(579.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCHO) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=O)CC1[CH]O1(15420)',
    structure = SMILES('[CH2]C(C=O)CC1[CH]O1'),
    E0 = (130.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.145788,0.0822817,-7.79122e-05,3.14857e-08,-1.13136e-12,15885.3,29.7853], Tmin=(100,'K'), Tmax=(818.38,'K')), NASAPolynomial(coeffs=[15.3621,0.0259179,-7.61227e-06,1.10761e-09,-6.56796e-14,12791.7,-44.2596], Tmin=(818.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.902,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(CJC(C)C=O)"""),
)

species(
    label = '[O]C=CCC1[CH]OC1(15421)',
    structure = SMILES('[O]C=CCC1[CH]OC1'),
    E0 = (64.2065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.616065,0.0547211,2.25613e-05,-8.73455e-08,4.46077e-11,7862.83,27.2668], Tmin=(100,'K'), Tmax=(883.103,'K')), NASAPolynomial(coeffs=[22.2016,0.0120796,1.349e-06,-6.40842e-10,5.00316e-14,1900.67,-86.3591], Tmin=(883.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.2065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Oxetane) + radical(C=COJ) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1[CH]CC(C=O)C1(15422)',
    structure = SMILES('[O]C1[CH]CC(C=O)C1'),
    E0 = (53.7784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53659,0.0359493,4.63379e-05,-8.44292e-08,3.43911e-11,6573.11,28.5583], Tmin=(100,'K'), Tmax=(970.419,'K')), NASAPolynomial(coeffs=[14.5882,0.0271499,-9.61593e-06,1.79411e-09,-1.3181e-13,1921.21,-44.9333], Tmin=(970.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.7784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1[CH]OOC=CC1(15385)',
    structure = SMILES('[CH2]C1[CH]OOC=CC1'),
    E0 = (324.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12871,0.00798727,0.000211507,-3.23995e-07,1.38562e-10,39148.1,28.7373], Tmin=(100,'K'), Tmax=(901.11,'K')), NASAPolynomial(coeffs=[40.7078,-0.0226015,2.08872e-05,-4.27309e-09,2.83014e-13,26123.9,-190.776], Tmin=(901.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCsJOOC) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(C=O)CC=CO(15423)',
    structure = SMILES('C=C(C=O)CC=CO'),
    E0 = (-250.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0687169,0.072891,-4.37126e-05,-2.46909e-09,7.67863e-12,-30016.7,29.4192], Tmin=(100,'K'), Tmax=(1011.62,'K')), NASAPolynomial(coeffs=[20.4015,0.0206592,-8.02695e-06,1.53531e-09,-1.1232e-13,-35571.6,-76.0326], Tmin=(1011.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'CC(C=O)CC=C=O(15424)',
    structure = SMILES('CC(C=O)CC=C=O'),
    E0 = (-237.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12941,0.0674545,-5.70318e-05,2.79367e-08,-6.02154e-12,-28501.4,27.5113], Tmin=(100,'K'), Tmax=(1057.57,'K')), NASAPolynomial(coeffs=[8.35368,0.0401309,-1.82782e-05,3.50784e-09,-2.4687e-13,-30029.5,-7.74614], Tmin=(1057.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH)"""),
)

species(
    label = '[CH2][C](C[O])CC=C[O](15425)',
    structure = SMILES('[CH2][C](C[O])CC=C[O]'),
    E0 = (283.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,360,370,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.387683,0.0714065,-5.76491e-05,2.46621e-08,-4.27743e-12,34183.4,34.0687], Tmin=(100,'K'), Tmax=(1372.06,'K')), NASAPolynomial(coeffs=[15.1091,0.0284896,-1.07314e-05,1.8658e-09,-1.23863e-13,30143.5,-41.6107], Tmin=(1372.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCOJ) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C[CH][O])C=O(15426)',
    structure = SMILES('[CH2]C([CH]C[CH][O])C=O'),
    E0 = (376.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0825221,0.0976874,-0.00015659,1.43873e-07,-5.1581e-11,45372,36.0076], Tmin=(100,'K'), Tmax=(843.085,'K')), NASAPolynomial(coeffs=[6.04532,0.0457311,-2.20446e-05,4.18662e-09,-2.86733e-13,45207.7,13.2464], Tmin=(843.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCJCC=O) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([CH]C=C[O])C[O](15427)',
    structure = SMILES('[CH2]C([CH]C=C[O])C[O]'),
    E0 = (271.619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.30322,0.0729623,-6.05131e-05,2.65308e-08,-4.6944e-12,32808.2,33.0698], Tmin=(100,'K'), Tmax=(1351.55,'K')), NASAPolynomial(coeffs=[15.57,0.0277793,-1.03673e-05,1.79574e-09,-1.19079e-13,28681.5,-45.1825], Tmin=(1351.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl) + radical(Allyl_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=C[O])CC[CH][O](15428)',
    structure = SMILES('[CH2]C(=C[O])CC[CH][O]'),
    E0 = (250.988,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,350,440,435,1725,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.126171,0.0843958,-8.5236e-05,4.5623e-08,-9.86809e-12,30327.1,31.6511], Tmin=(100,'K'), Tmax=(1112.74,'K')), NASAPolynomial(coeffs=[14.9921,0.0309562,-1.31973e-05,2.46242e-09,-1.71049e-13,27018.7,-41.6561], Tmin=(1112.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CCOJ) + radical(CCsJOH) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH]O)[CH]C=C[O](15429)',
    structure = SMILES('[CH2]C([CH]O)[CH]C=C[O]'),
    E0 = (226.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.505767,0.0866094,-8.73821e-05,4.47896e-08,-8.88251e-12,27379.8,34.658], Tmin=(100,'K'), Tmax=(1307.33,'K')), NASAPolynomial(coeffs=[20.3755,0.0189522,-5.43086e-06,7.94566e-10,-4.78081e-14,22242,-70.4464], Tmin=(1307.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOH) + radical(Allyl_S) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C[O])C[CH]C[O](15430)',
    structure = SMILES('[CH2]C(=C[O])C[CH]C[O]'),
    E0 = (270.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,350,440,435,1725,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.103945,0.0760876,-6.46943e-05,2.81469e-08,-4.88236e-12,32692.9,34.6826], Tmin=(100,'K'), Tmax=(1385.23,'K')), NASAPolynomial(coeffs=[17.6473,0.0254286,-9.83752e-06,1.74584e-09,-1.17557e-13,27832.6,-55.6708], Tmin=(1385.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CCJCO) + radical(C=COJ) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)CC[CH][O](15431)',
    structure = SMILES('[CH2]C([C]=O)CC[CH][O]'),
    E0 = (334.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,1855,455,950,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0523636,0.101542,-0.000166998,1.54142e-07,-5.50706e-11,40419.2,35.152], Tmin=(100,'K'), Tmax=(849.595,'K')), NASAPolynomial(coeffs=[6.28129,0.0455798,-2.20403e-05,4.17816e-09,-2.85303e-13,40286.5,11.1807], Tmin=(849.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOH) + radical(CC(C)CJ=O) + radical(CCOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C[O])C[C]=C[O](15432)',
    structure = SMILES('[CH2]C(C[O])C[C]=C[O]'),
    E0 = (368.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469884,0.0781186,-7.81155e-05,4.37531e-08,-1.00848e-11,44428.8,34.0521], Tmin=(100,'K'), Tmax=(1040.12,'K')), NASAPolynomial(coeffs=[12.0726,0.0334984,-1.37671e-05,2.50914e-09,-1.71583e-13,42015.1,-22.3806], Tmin=(1040.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([C]=O)C[CH]C[O](15433)',
    structure = SMILES('[CH2]C([C]=O)C[CH]C[O]'),
    E0 = (354.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,1855,455,950,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.401687,0.0877031,-0.000127601,1.13061e-07,-4.04462e-11,42764.4,36.4701], Tmin=(100,'K'), Tmax=(815.802,'K')), NASAPolynomial(coeffs=[6.22038,0.0444972,-2.11746e-05,4.03901e-09,-2.78989e-13,42303.4,12.5761], Tmin=(815.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CJC(C)C=O) + radical(CC(C)CJ=O) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([CH]O)C[C]=C[O](15434)',
    structure = SMILES('[CH2]C([CH]O)C[C]=C[O]'),
    E0 = (322.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.360671,0.0920481,-0.000106043,6.3379e-08,-1.48198e-11,39001.2,35.715], Tmin=(100,'K'), Tmax=(1052.88,'K')), NASAPolynomial(coeffs=[17.7461,0.0232585,-8.04002e-06,1.32511e-09,-8.53803e-14,35188.3,-52.5725], Tmin=(1052.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOH) + radical(Cds_S) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C[O])C[CH][C]=O(15435)',
    structure = SMILES('[CH2]C(C[O])C[CH][C]=O'),
    E0 = (318.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,1855,455,950,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541677,0.0820892,-0.000108488,8.88047e-08,-2.96133e-11,38376.3,34.1849], Tmin=(100,'K'), Tmax=(850.419,'K')), NASAPolynomial(coeffs=[7.05049,0.0410484,-1.77091e-05,3.22389e-09,-2.1677e-13,37646.3,6.0549], Tmin=(850.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCCJ=O) + radical(CCOJ) + radical(CCJCHO)"""),
)

species(
    label = '[CH2]C([CH]O)C[CH][C]=O(15436)',
    structure = SMILES('[CH2]C([CH]O)C[CH][C]=O'),
    E0 = (272.685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0865514,0.0935426,-0.000127288,9.58025e-08,-2.85473e-11,32940,35.1277], Tmin=(100,'K'), Tmax=(896.674,'K')), NASAPolynomial(coeffs=[12.3522,0.0314444,-1.23499e-05,2.12693e-09,-1.37797e-13,30975,-22.044], Tmin=(896.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCJCHO) + radical(CCsJOH) + radical(CCCJ=O)"""),
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
    label = '[CH2]CCC=C[O](6414)',
    structure = SMILES('[CH2]CCC=C[O]'),
    E0 = (97.6143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,364.629,364.63,364.64],'cm^-1')),
        HinderedRotor(inertia=(0.00126791,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194436,'amu*angstrom^2'), symmetry=1, barrier=(18.344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194428,'amu*angstrom^2'), symmetry=1, barrier=(18.344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3713.27,'J/mol'), sigma=(6.31777,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=580.00 K, Pc=33.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3357,0.0455502,1.93912e-06,-3.8633e-08,1.90622e-11,11848,24.941], Tmin=(100,'K'), Tmax=(963.177,'K')), NASAPolynomial(coeffs=[15.2034,0.0194655,-6.50506e-06,1.17361e-09,-8.50461e-14,7715.12,-49.029], Tmin=(963.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.6143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CCCC=C[O](12912)',
    structure = SMILES('[O]C=CCCC=C[O]'),
    E0 = (-70.9601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266561,0.0614179,4.07998e-06,-6.33251e-08,3.25859e-11,-8381.13,30.2869], Tmin=(100,'K'), Tmax=(940.103,'K')), NASAPolynomial(coeffs=[23.7137,0.0137017,-2.83136e-06,4.673e-10,-3.91771e-14,-15089.6,-93.6163], Tmin=(940.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.9601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CC[CH]CC=O(15437)',
    structure = SMILES('[O]C=CC[CH]CC=O'),
    E0 = (-14.3096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.454382,0.0675237,-4.19804e-05,4.79243e-09,3.27901e-12,-1584.39,33.0928], Tmin=(100,'K'), Tmax=(1049.38,'K')), NASAPolynomial(coeffs=[16.2063,0.0264354,-1.03417e-05,1.90473e-09,-1.33579e-13,-5933.97,-48.6329], Tmin=(1049.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.3096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCJCC=O)"""),
)

species(
    label = 'O=CC1CC=COC1(15438)',
    structure = SMILES('O=CC1CC=COC1'),
    E0 = (-266.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31433,0.0386327,5.05859e-05,-9.93037e-08,4.30138e-11,-31964.3,21.278], Tmin=(100,'K'), Tmax=(930.681,'K')), NASAPolynomial(coeffs=[17.5616,0.020983,-5.0672e-06,8.0416e-10,-5.95051e-14,-37248.3,-68.0794], Tmin=(930.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran)"""),
)

species(
    label = '[CH2]C=COCC=C[O](12913)',
    structure = SMILES('[CH2]C=COCC=C[O]'),
    E0 = (-27.8045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0453551,0.0661676,1.76124e-06,-6.89979e-08,3.66818e-11,-3177.4,29.6622], Tmin=(100,'K'), Tmax=(930.055,'K')), NASAPolynomial(coeffs=[26.8323,0.00898402,-2.20484e-07,-4.88871e-11,-3.47383e-15,-10703.3,-111.64], Tmin=(930.055,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.8045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(=CO)CC=C[O](15439)',
    structure = SMILES('[CH2]C(=CO)CC=C[O]'),
    E0 = (-74.9287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.09668,'amu*angstrom^2'), symmetry=1, barrier=(25.2149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09682,'amu*angstrom^2'), symmetry=1, barrier=(25.218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09697,'amu*angstrom^2'), symmetry=1, barrier=(25.2214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09701,'amu*angstrom^2'), symmetry=1, barrier=(25.2223,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254689,0.0699241,-2.94804e-06,-7.03467e-08,3.90964e-11,-8836.64,30.2619], Tmin=(100,'K'), Tmax=(917.918,'K')), NASAPolynomial(coeffs=[28.9985,0.00452591,2.47746e-06,-6.10335e-10,3.70053e-14,-16822.3,-122.608], Tmin=(917.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.9287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]OC([O])[CH]C1(15440)',
    structure = SMILES('[CH2]C1[CH]OC([O])[CH]C1'),
    E0 = (333.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24279,0.0481847,1.5638e-05,-5.96021e-08,2.93388e-11,40273.6,27.7291], Tmin=(100,'K'), Tmax=(887.283,'K')), NASAPolynomial(coeffs=[13.6446,0.0267847,-6.52334e-06,8.82489e-10,-5.38204e-14,36714.4,-38.2739], Tmin=(887.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxane) + radical(Isobutyl) + radical(CCJCO) + radical(CCOJ) + radical(CCsJOCs)"""),
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
    label = '[O]C=CCC=C[O](14266)',
    structure = SMILES('[O]C=CCC=C[O]'),
    E0 = (-45.9103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15021,'amu*angstrom^2'), symmetry=1, barrier=(26.4457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15132,'amu*angstrom^2'), symmetry=1, barrier=(26.4712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.988195,0.0456095,2.17824e-05,-7.67715e-08,3.71802e-11,-5394.22,25.4798], Tmin=(100,'K'), Tmax=(928.446,'K')), NASAPolynomial(coeffs=[22.5849,0.0055759,8.15949e-07,-2.19695e-10,8.27632e-15,-11689.3,-89.4131], Tmin=(928.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.9103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[CH]=CCC([CH2])C=O(15441)',
    structure = SMILES('[CH]=CCC([CH2])C=O'),
    E0 = (304.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180.552],'cm^-1')),
        HinderedRotor(inertia=(0.188624,'amu*angstrom^2'), symmetry=1, barrier=(4.33683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574095,'amu*angstrom^2'), symmetry=1, barrier=(13.1996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.029464,'amu*angstrom^2'), symmetry=1, barrier=(10.7275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.461146,'amu*angstrom^2'), symmetry=1, barrier=(10.6054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01634,0.0695769,-7.05646e-05,4.17497e-08,-1.05047e-11,36782.9,27.939], Tmin=(100,'K'), Tmax=(940.027,'K')), NASAPolynomial(coeffs=[9.05077,0.0353885,-1.60094e-05,3.05873e-09,-2.14672e-13,35272.4,-10.3255], Tmin=(940.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(C=O)CC=C[O](15442)',
    structure = SMILES('[CH]C(C=O)CC=C[O]'),
    E0 = (228.243,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209402,0.0738367,-6.04384e-05,2.05981e-08,-1.30172e-12,27596,32.3408], Tmin=(100,'K'), Tmax=(1046.18,'K')), NASAPolynomial(coeffs=[17.7442,0.0221862,-8.45231e-06,1.53431e-09,-1.06895e-13,23084.7,-57.0719], Tmin=(1046.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC(C=O)C1[O](15443)',
    structure = SMILES('[CH2]C1CC(C=O)C1[O]'),
    E0 = (131.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04708,0.0512739,3.79726e-06,-4.11919e-08,1.93954e-11,15956.4,28.4353], Tmin=(100,'K'), Tmax=(980.092,'K')), NASAPolynomial(coeffs=[14.5679,0.028272,-1.02494e-05,1.86324e-09,-1.32253e-13,11760.5,-44.4076], Tmin=(980.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(Isobutyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(C=O)C=CC=O(15444)',
    structure = SMILES('[CH2]C(C=O)C=CC=O'),
    E0 = (-42.974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.247266,'amu*angstrom^2'), symmetry=1, barrier=(7.14536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.496404,'amu*angstrom^2'), symmetry=1, barrier=(11.4133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.495784,'amu*angstrom^2'), symmetry=1, barrier=(11.399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.495916,'amu*angstrom^2'), symmetry=1, barrier=(11.4021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.556035,0.0836971,-0.000119792,1.0674e-07,-3.96618e-11,-5052.21,29.3366], Tmin=(100,'K'), Tmax=(748.919,'K')), NASAPolynomial(coeffs=[6.39813,0.0437668,-2.23363e-05,4.42705e-09,-3.13959e-13,-5682.51,4.47494], Tmin=(748.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([CH]CC=O)C=O(15321)',
    structure = SMILES('[CH2]C([CH]CC=O)C=O'),
    E0 = (47.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,240.069,1153.59],'cm^-1')),
        HinderedRotor(inertia=(0.00292514,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132605,'amu*angstrom^2'), symmetry=1, barrier=(5.4247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13261,'amu*angstrom^2'), symmetry=1, barrier=(5.42344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132606,'amu*angstrom^2'), symmetry=1, barrier=(5.42364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132673,'amu*angstrom^2'), symmetry=1, barrier=(5.42346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23574,0.0714866,-5.04281e-05,-3.01349e-08,5.38959e-11,5764.28,30.9418], Tmin=(100,'K'), Tmax=(499.991,'K')), NASAPolynomial(coeffs=[7.13952,0.0428337,-2.02038e-05,3.88112e-09,-2.71044e-13,4941.69,4.22949], Tmin=(499.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=O)CC[C]=O(12901)',
    structure = SMILES('[CH2]C(C=O)CC[C]=O'),
    E0 = (7.2475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,289.716,2185.41],'cm^-1')),
        HinderedRotor(inertia=(0.15355,'amu*angstrom^2'), symmetry=1, barrier=(9.14581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15355,'amu*angstrom^2'), symmetry=1, barrier=(9.14581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15355,'amu*angstrom^2'), symmetry=1, barrier=(9.14581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15355,'amu*angstrom^2'), symmetry=1, barrier=(9.14581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15355,'amu*angstrom^2'), symmetry=1, barrier=(9.14581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.212663,0.0914172,-0.000134091,1.15713e-07,-4.00722e-11,1000.24,32.783], Tmin=(100,'K'), Tmax=(825.564,'K')), NASAPolynomial(coeffs=[7.76292,0.0415749,-1.94377e-05,3.67187e-09,-2.51965e-13,205.462,0.541408], Tmin=(825.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.2475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(=C[O])CCC=O(15445)',
    structure = SMILES('[CH2]C(=C[O])CCC=O'),
    E0 = (-77.9873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.145592,0.074461,-5.51293e-05,1.4978e-08,4.19165e-13,-9232.05,30.4974], Tmin=(100,'K'), Tmax=(1053.21,'K')), NASAPolynomial(coeffs=[17.5491,0.0255251,-9.87496e-06,1.8033e-09,-1.25794e-13,-13849.8,-58.8854], Tmin=(1053.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.9873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)CCC=O(15446)',
    structure = SMILES('[CH2]C([C]=O)CCC=O'),
    E0 = (5.97554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,201.409,3416.18],'cm^-1')),
        HinderedRotor(inertia=(0.00419881,'amu*angstrom^2'), symmetry=1, barrier=(0.12465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16353,'amu*angstrom^2'), symmetry=1, barrier=(34.1855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31602,'amu*angstrom^2'), symmetry=1, barrier=(9.11647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.313024,'amu*angstrom^2'), symmetry=1, barrier=(9.10654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30844,'amu*angstrom^2'), symmetry=1, barrier=(9.13129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.37569,0.0867271,-0.000119582,1.00795e-07,-3.50418e-11,842.471,32.5369], Tmin=(100,'K'), Tmax=(790.406,'K')), NASAPolynomial(coeffs=[7.54163,0.0423559,-1.99915e-05,3.81992e-09,-2.65006e-13,-37.1028,1.25271], Tmin=(790.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.97554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CC(C)CJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C1[CH]OC(C=O)C1(15344)',
    structure = SMILES('[CH2]C1[CH]OC(C=O)C1'),
    E0 = (37.5208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29008,0.0387271,5.6768e-05,-1.14604e-07,5.22576e-11,4630.02,26.1339], Tmin=(100,'K'), Tmax=(889.124,'K')), NASAPolynomial(coeffs=[18.8825,0.0160814,-3.44497e-07,-3.11816e-10,2.61648e-14,-731.587,-69.2305], Tmin=(889.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(Isobutyl)"""),
)

species(
    label = 'O=CC1C[CH][CH]OC1(15447)',
    structure = SMILES('O=CC1C[CH][CH]OC1'),
    E0 = (23.8554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22129,0.047045,2.10977e-05,-6.65604e-08,3.18986e-11,2982.4,24.8701], Tmin=(100,'K'), Tmax=(896.957,'K')), NASAPolynomial(coeffs=[14.6469,0.0251966,-5.95213e-06,8.05919e-10,-5.04277e-14,-955.59,-46.967], Tmin=(896.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.8554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Oxane) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'C=C(C=O)CCC=O(15448)',
    structure = SMILES('C=C(C=O)CCC=O'),
    E0 = (-253.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28466,0.0673858,-5.89068e-05,3.25019e-08,-8.57907e-12,-30446.8,26.7579], Tmin=(100,'K'), Tmax=(841.174,'K')), NASAPolynomial(coeffs=[5.40638,0.0477858,-2.39555e-05,4.80127e-09,-3.46297e-13,-31140.2,7.58596], Tmin=(841.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'CC(C=O)C=CC=O(15449)',
    structure = SMILES('CC(C=O)C=CC=O'),
    E0 = (-253.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18825,0.0690942,-6.53758e-05,4.04436e-08,-1.18238e-11,-30392.3,26.7484], Tmin=(100,'K'), Tmax=(777.075,'K')), NASAPolynomial(coeffs=[5.37442,0.0475456,-2.37797e-05,4.75714e-09,-3.4262e-13,-31042.9,7.60837], Tmin=(777.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = '[CH2]C(C=O)C([CH2])C=O(12917)',
    structure = SMILES('[CH2]C(C=O)C([CH2])C=O'),
    E0 = (52.0371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.200371,'amu*angstrom^2'), symmetry=1, barrier=(4.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20059,'amu*angstrom^2'), symmetry=1, barrier=(4.61195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199605,'amu*angstrom^2'), symmetry=1, barrier=(4.5893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20081,'amu*angstrom^2'), symmetry=1, barrier=(4.61701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.702052,'amu*angstrom^2'), symmetry=1, barrier=(16.1416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3978.28,'J/mol'), sigma=(6.59354,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=621.40 K, Pc=31.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00699609,0.096428,-0.00014472,1.2514e-07,-4.32009e-11,6394.12,32.2225], Tmin=(100,'K'), Tmax=(825.648,'K')), NASAPolynomial(coeffs=[8.47563,0.0414858,-1.96242e-05,3.72011e-09,-2.55557e-13,5469.97,-4.13919], Tmin=(825.648,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.0371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CJC(C)C=O)"""),
)

species(
    label = 'O=CC1CC(C=O)C1(12920)',
    structure = SMILES('O=CC1CC(C=O)C1'),
    E0 = (-210.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26791,0.0487367,-8.62281e-07,-2.68708e-08,1.16733e-11,-25246.3,25.181], Tmin=(100,'K'), Tmax=(1072.14,'K')), NASAPolynomial(coeffs=[11.9012,0.0336008,-1.40128e-05,2.65096e-09,-1.8769e-13,-28936.5,-33.4353], Tmin=(1072.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH]CC([CH2])C=O(12942)',
    structure = SMILES('[CH]CC([CH2])C=O'),
    E0 = (402.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180.163,697.062,2610,2611.18],'cm^-1')),
        HinderedRotor(inertia=(0.549407,'amu*angstrom^2'), symmetry=1, barrier=(12.6319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00717067,'amu*angstrom^2'), symmetry=1, barrier=(2.45816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87997,'amu*angstrom^2'), symmetry=1, barrier=(43.2242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549492,'amu*angstrom^2'), symmetry=1, barrier=(12.6339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3289,0.0626831,-7.1191e-05,4.82373e-08,-1.3803e-11,48451.9,24.6975], Tmin=(100,'K'), Tmax=(836.678,'K')), NASAPolynomial(coeffs=[8.00286,0.0307773,-1.39926e-05,2.66328e-09,-1.85943e-13,47335.1,-6.31075], Tmin=(836.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CJC(C)C=O)"""),
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
    E0 = (-9.46149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (80.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (114.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (112.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (113.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (187.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (60.0687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (74.4419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (118.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (251.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (109.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (108.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (237.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (62.1781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (104.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (19.0013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (78.8209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (92.5448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (180.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (325.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (278.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (410.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (407.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (440.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (361.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (386.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (204.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (176.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (53.7784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (324.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (15.5118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (15.5118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (305.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (398.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (360.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (339.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (265.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (309.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (359.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (407.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (379.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (362.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (343.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (297.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (702.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (150.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (113.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-1.26085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (285.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (68.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (333.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (369.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (711.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (440.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (131.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (168.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (172.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (166.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (191.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (108.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (146.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (78.6208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (79.5071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (68.7856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (211.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (-1.17717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (469.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['C=CC=O(5269)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[O]C=CCC1CC1[O](15400)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.93521e+09,'s^-1'), n=0.743095, Ea=(89.4945,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 85.9 to 89.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[O][CH]C1CC(C=O)C1(15401)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(177207,'s^-1'), n=1.88643, Ea=(124.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 121.6 to 124.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[CH2]C1CC=COC1[O](15402)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(121.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;multiplebond_intra;radadd_intra_O] for rate rule [R7_SMSS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C=C(C=O)CC=C[O](15403)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(72.1434,'m^3/(mol*s)'), n=1.66666, Ea=(10.8177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;HJ] for rate rule [Cds-COCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', '[CH2]C(C=O)CC=C=O(15404)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C=C[O](5266)', 'C=CC=O(5269)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0123754,'m^3/(mol*s)'), n=2.41, Ea=(51.1145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeH_Cds;CsJ-CdHH] for rate rule [Cds-COH_Cds;CsJ-CdHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=O(373)', 'C=CCC=C[O](12633)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CJ] for rate rule [Cds-CsH_Cds-HH;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC(=C[O])CC=C[O](15405)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C=O)CC=[C]O(15406)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['CC([CH]C=C[O])C=O(15407)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(25000,'s^-1'), n=2.28, Ea=(119.244,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 85 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['CC([C]=O)CC=C[O](15408)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_2H;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=O)C[C]=CO(15409)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CC(C=O)C[C]=C[O](15410)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[CH2]C(C=O)C=C[CH]O(15411)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 289 used for R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['CC(C=O)C[CH][C]=O(15412)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.169e+07,'s^-1'), n=1.0878, Ea=(28.4628,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSR;C_rad_out_2H;XH_out] for rate rule [R5H_SSSD;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[CH2]C(=C[O])CC=CO(15413)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;XH_out] for rate rule [R5H_SMSS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([C]=O)CC=CO(15414)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(16140,'s^-1'), n=1.92259, Ea=(84.7802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C=C[O](5266)', '[CH2]C=C[O](5266)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.13324e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=O(373)', '[CH2][CH]CC=C[O](6404)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.7839e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2]C(=C[O])CC=C[O](15415)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C[O](602)', '[CH2]C([CH2])C=O(6124)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.00436e+08,'m^3/(mol*s)'), n=-0.446058, Ea=(0.74957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H2/Cs] + [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/Cs]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2]C([CH]C=C[O])C=O(15416)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.66284e+07,'m^3/(mol*s)'), n=0.108445, Ea=(5.74998,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH2]C(C=O)C[C]=C[O](15417)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH2]C([C]=O)CC=C[O](15418)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH2]C(C=O)C[CH][C]=O(15419)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[CH2]C(C=O)CC1[CH]O1(15420)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[O]C=CCC1[CH]OC1(15421)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[O]C1[CH]CC(C=O)C1(15422)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.89094e+07,'s^-1'), n=0.979167, Ea=(63.2399,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_CsCs_RH_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 60.1 to 63.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[CH2]C1[CH]OOC=CC1(15385)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(333.666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 325.4 to 333.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['C=C(C=O)CC=CO(15423)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['CC(C=O)CC=C=O(15424)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][C](C[O])CC=C[O](15425)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH]C[CH][O])C=O(15426)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([CH]C=C[O])C[O](15427)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(=C[O])CC[CH][O](15428)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH]O)[CH]C=C[O](15429)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(=C[O])C[CH]C[O](15430)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([C]=O)CC[CH][O](15431)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(C[O])C[C]=C[O](15432)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([C]=O)C[CH]C[O](15433)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([CH]O)C[C]=C[O](15434)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C(C[O])C[CH][C]=O(15435)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([CH]O)C[CH][C]=O(15436)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[C-]#[O+](374)', '[CH2]CCC=C[O](6414)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[O]C=CCCC=C[O](12912)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[O]C=CC[CH]CC=O(15437)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['O=CC1CC=COC1(15438)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Ypri_rad_out] + [R6_SSSDS;C_rad_out_single;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C=COCC=C[O](12913)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C(=CO)CC=C[O](15439)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C1[CH]OC([O])[CH]C1(15440)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction52',
    reactants = ['CH2(T)(28)', '[O]C=CCC=C[O](14266)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction53',
    reactants = ['O(T)(63)', '[CH]=CCC([CH2])C=O(15441)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction54',
    reactants = ['H(8)', '[CH]C(C=O)CC=C[O](15442)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[CH2]C1CC(C=O)C1[O](15443)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(413769,'s^-1'), n=1.87624, Ea=(141.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 137.3 to 141.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction56',
    reactants = ['H(8)', '[CH2]C(C=O)C=CC=O(15444)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(4.76955,'m^3/(mol*s)'), n=1.94497, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-OneDeH;HJ] for rate rule [Cds-CsH_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C([CH]CC=O)C=O(15321)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2]C(C=O)CC[C]=O(12901)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[CH2]C(=C[O])CCC=O(15445)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(3.07201e+08,'s^-1'), n=1.25033, Ea=(201.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/OneDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C([C]=O)CCC=O(15446)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;CO_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['[CH2]C1[CH]OC(C=O)C1(15344)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHCO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['O=CC1C[CH][CH]OC1(15447)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(88.0823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['C=C(C=O)CCC=O(15448)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['CC(C=O)C=CC=O(15449)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=O)C([CH2])C=O(12917)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C(C=O)CC=C[O](12914)'],
    products = ['O=CC1CC(C=O)C1(12920)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH]=O(373)', '[CH]CC([CH2])C=O(12942)'],
    products = ['[CH2]C(C=O)CC=C[O](12914)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '3676',
    isomers = [
        '[CH2]C(C=O)CC=C[O](12914)',
    ],
    reactants = [
        ('C=CC=O(5269)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3676',
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

