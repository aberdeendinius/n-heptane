species(
    label = 'C=[C][CH]C([O])C=C(18289)',
    structure = SMILES('C=[C][CH]C([O])C=C'),
    E0 = (429.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,546.233,546.346,546.388,546.508],'cm^-1')),
        HinderedRotor(inertia=(0.12446,'amu*angstrom^2'), symmetry=1, barrier=(26.3541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12436,'amu*angstrom^2'), symmetry=1, barrier=(26.3462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124248,'amu*angstrom^2'), symmetry=1, barrier=(26.341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20991,0.0447833,1.47304e-05,-5.50627e-08,2.49446e-11,51740.9,29.403], Tmin=(100,'K'), Tmax=(978.524,'K')), NASAPolynomial(coeffs=[17.2383,0.019111,-7e-06,1.35861e-09,-1.02782e-13,46696.3,-57.3247], Tmin=(978.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ) + radical(C=CCJC(O)C=C)"""),
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
    label = '[CH2]C1OC1[CH][C]=C(20294)',
    structure = SMILES('[CH2]C1OC1[CH][C]=C'),
    E0 = (506.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.154879,0.0742211,-7.15048e-05,3.57327e-08,-6.75914e-12,61108.1,26.1465], Tmin=(100,'K'), Tmax=(1508.08,'K')), NASAPolynomial(coeffs=[16.932,0.0170506,-2.85449e-06,1.74733e-10,-8.34001e-16,57301.9,-58.8404], Tmin=(1508.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(CJCO) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1C([O])C1[C]=C(20134)',
    structure = SMILES('[CH2]C1C([O])C1[C]=C'),
    E0 = (582.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36099,0.0440111,1.03138e-05,-5.01647e-08,2.41721e-11,70144.2,26.1028], Tmin=(100,'K'), Tmax=(939.126,'K')), NASAPolynomial(coeffs=[15.6383,0.0184757,-5.24403e-06,8.76705e-10,-6.28747e-14,65907,-50.1625], Tmin=(939.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S) + radical(Isobutyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1C(=C)[CH]C1[O](20100)',
    structure = SMILES('[CH2]C1C(=C)[CH]C1[O]'),
    E0 = (448.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.5398,-0.0103592,0.000121655,-1.19437e-07,2.76053e-11,53588.7,-18.5974], Tmin=(100,'K'), Tmax=(1723.24,'K')), NASAPolynomial(coeffs=[76.64,0.0222368,-6.86475e-05,1.68301e-08,-1.25196e-12,3186.28,-453.608], Tmin=(1723.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C=CCJCO) + radical(Isobutyl) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C=C([O])C=C=C(20295)',
    structure = SMILES('[CH2]C=C([O])C=C=C'),
    E0 = (240.442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.4319,'amu*angstrom^2'), symmetry=1, barrier=(32.9222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43197,'amu*angstrom^2'), symmetry=1, barrier=(32.9239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980007,0.0582019,-3.7442e-05,8.20562e-10,5.95249e-12,29034.7,23.802], Tmin=(100,'K'), Tmax=(942.213,'K')), NASAPolynomial(coeffs=[15.0929,0.0182996,-5.77614e-06,9.56903e-10,-6.47259e-14,25487,-48.1584], Tmin=(942.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C=CC([O])C=C(20057)',
    structure = SMILES('[CH]=C=CC([O])C=C'),
    E0 = (436.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.970908,'amu*angstrom^2'), symmetry=1, barrier=(22.3231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.9706,'amu*angstrom^2'), symmetry=1, barrier=(22.316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0173,0.0582997,-4.56264e-05,1.57995e-08,-1.35038e-12,52592.3,27.902], Tmin=(100,'K'), Tmax=(1075.48,'K')), NASAPolynomial(coeffs=[13.8762,0.0208721,-7.92724e-06,1.42018e-09,-9.74795e-14,49225.1,-37.8661], Tmin=(1075.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C=CC=C=C(19683)',
    structure = SMILES('[CH2]C=CC=C=C'),
    E0 = (316.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.57119,'amu*angstrom^2'), symmetry=1, barrier=(36.1247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57201,'amu*angstrom^2'), symmetry=1, barrier=(36.1435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67678,0.0391981,9.13911e-06,-4.16281e-08,1.94046e-11,38164.2,20.5572], Tmin=(100,'K'), Tmax=(956.834,'K')), NASAPolynomial(coeffs=[13.1801,0.0202038,-6.69503e-06,1.18327e-09,-8.42569e-14,34631,-41.3917], Tmin=(956.834,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=C(64)',
    structure = SMILES('[CH]=C'),
    E0 = (289.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,826.012,826.012,3240.27],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90671,-0.00406241,3.8678e-05,-4.62976e-08,1.729e-11,34797.2,6.09789], Tmin=(100,'K'), Tmax=(931.962,'K')), NASAPolynomial(coeffs=[5.44797,0.00498356,-1.08821e-06,1.79837e-10,-1.45096e-14,33829.8,-4.87808], Tmin=(931.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C=C[O](18052)',
    structure = SMILES('C=[C]C=C[O]'),
    E0 = (225.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61747,'amu*angstrom^2'), symmetry=1, barrier=(37.1889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98082,0.0322081,7.54003e-06,-4.54908e-08,2.38758e-11,27216,16.1105], Tmin=(100,'K'), Tmax=(899.941,'K')), NASAPolynomial(coeffs=[16.1068,0.00249527,1.93911e-06,-5.05291e-10,3.47522e-14,23334.2,-57.9907], Tmin=(899.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CC(O)=C[CH2](20296)',
    structure = SMILES('[CH2][C]=CC(O)=C[CH2]'),
    E0 = (281.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835691,0.0580437,-2.16502e-05,-2.25839e-08,1.59892e-11,34033,25.798], Tmin=(100,'K'), Tmax=(919.664,'K')), NASAPolynomial(coeffs=[17.1659,0.0172447,-4.40847e-06,6.57182e-10,-4.40493e-14,29751.1,-58.568], Tmin=(919.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=C([O])C[C]=C(18288)',
    structure = SMILES('[CH2]C=C([O])C[C]=C'),
    E0 = (365.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,378.552,379.222,379.399],'cm^-1')),
        HinderedRotor(inertia=(0.00116957,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143982,'amu*angstrom^2'), symmetry=1, barrier=(14.7702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144176,'amu*angstrom^2'), symmetry=1, barrier=(14.7838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86388,0.0611858,-4.54149e-05,1.25769e-08,3.89637e-13,44056.5,28.1986], Tmin=(100,'K'), Tmax=(1025.45,'K')), NASAPolynomial(coeffs=[14.5012,0.0218156,-8.04851e-06,1.43163e-09,-9.84796e-14,40532.7,-41.4805], Tmin=(1025.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C[CH]C([O])C=C(20059)',
    structure = SMILES('[CH]=C[CH]C([O])C=C'),
    E0 = (438.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,500.923,500.932,501.076],'cm^-1')),
        HinderedRotor(inertia=(0.155645,'amu*angstrom^2'), symmetry=1, barrier=(27.7149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15559,'amu*angstrom^2'), symmetry=1, barrier=(27.714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155656,'amu*angstrom^2'), symmetry=1, barrier=(27.7144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16968,0.0437024,2.28644e-05,-6.69423e-08,2.99871e-11,52857.2,29.4766], Tmin=(100,'K'), Tmax=(967.882,'K')), NASAPolynomial(coeffs=[18.6374,0.0168288,-5.71705e-06,1.11754e-09,-8.71367e-14,47353.3,-65.1899], Tmin=(967.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C][CH]C(O)[C]=C(20297)',
    structure = SMILES('C=[C][CH]C(O)[C]=C'),
    E0 = (436.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1277.5,1000,1380,1390,370,380,2900,435,705.376,705.387,705.395],'cm^-1')),
        HinderedRotor(inertia=(0.480567,'amu*angstrom^2'), symmetry=1, barrier=(20.8155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905336,'amu*angstrom^2'), symmetry=1, barrier=(20.8155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0589514,'amu*angstrom^2'), symmetry=1, barrier=(20.8154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.058954,'amu*angstrom^2'), symmetry=1, barrier=(20.8154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00748,0.0534506,-1.61976e-05,-2.04947e-08,1.22809e-11,52644,30.7954], Tmin=(100,'K'), Tmax=(999.958,'K')), NASAPolynomial(coeffs=[16.1987,0.020262,-7.78247e-06,1.47611e-09,-1.07528e-13,48227.1,-49.3871], Tmin=(999.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC([O])[C]=C(18290)',
    structure = SMILES('C=[C]CC([O])[C]=C'),
    E0 = (596.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,337.232,337.287,337.292,4000],'cm^-1')),
        HinderedRotor(inertia=(0.237396,'amu*angstrom^2'), symmetry=1, barrier=(19.1604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237331,'amu*angstrom^2'), symmetry=1, barrier=(19.1606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23734,'amu*angstrom^2'), symmetry=1, barrier=(19.1603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13791,0.0609972,-5.09129e-05,2.22672e-08,-3.99731e-12,71899.5,29.3365], Tmin=(100,'K'), Tmax=(1307.01,'K')), NASAPolynomial(coeffs=[12.368,0.0266285,-1.14694e-05,2.14824e-09,-1.49045e-13,68963.9,-27.849], Tmin=(1307.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=CC([O])=C[CH2](20298)',
    structure = SMILES('[CH2]C=CC([O])=C[CH2]'),
    E0 = (181.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22076,0.0496472,-5.05699e-06,-3.37748e-08,1.84451e-11,21987.4,25.4427], Tmin=(100,'K'), Tmax=(928.425,'K')), NASAPolynomial(coeffs=[14.6638,0.0210239,-6.14143e-06,9.89231e-10,-6.72696e-14,18228.6,-45.2135], Tmin=(928.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=[C]CC([O])C=C(18292)',
    structure = SMILES('[CH]=[C]CC([O])C=C'),
    E0 = (606.192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,310.737,311.378,312.424],'cm^-1')),
        HinderedRotor(inertia=(0.00172324,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22347,'amu*angstrom^2'), symmetry=1, barrier=(15.3409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223474,'amu*angstrom^2'), symmetry=1, barrier=(15.3497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.5,'J/mol'), sigma=(6.3848,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.16 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.845992,0.0627885,-5.23916e-05,2.22753e-08,-3.79546e-12,73027,30.3196], Tmin=(100,'K'), Tmax=(1399.64,'K')), NASAPolynomial(coeffs=[14.9752,0.0224084,-9.11566e-06,1.66213e-09,-1.13562e-13,69071.9,-42.5962], Tmin=(1399.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C([O])[CH]C=C(20299)',
    structure = SMILES('C=[C]C([O])[CH]C=C'),
    E0 = (429.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,527.406,527.647,528.058,529.206],'cm^-1')),
        HinderedRotor(inertia=(0.127981,'amu*angstrom^2'), symmetry=1, barrier=(24.9895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126621,'amu*angstrom^2'), symmetry=1, barrier=(25.0487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.598911,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2099,0.0447834,1.47302e-05,-5.50624e-08,2.49444e-11,51740.9,29.403], Tmin=(100,'K'), Tmax=(978.525,'K')), NASAPolynomial(coeffs=[17.2383,0.019111,-6.99998e-06,1.35861e-09,-1.02782e-13,46696.3,-57.3248], Tmin=(978.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(O)[CH][C]=C(20300)',
    structure = SMILES('[CH]=CC(O)[CH][C]=C'),
    E0 = (445.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,388.895,389.455],'cm^-1')),
        HinderedRotor(inertia=(0.245568,'amu*angstrom^2'), symmetry=1, barrier=(26.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111417,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244999,'amu*angstrom^2'), symmetry=1, barrier=(26.2217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244898,'amu*angstrom^2'), symmetry=1, barrier=(26.2312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972877,0.0523075,-7.87135e-06,-3.2579e-08,1.73874e-11,53760.1,30.8486], Tmin=(100,'K'), Tmax=(980.407,'K')), NASAPolynomial(coeffs=[17.547,0.0180641,-6.5473e-06,1.24618e-09,-9.27981e-14,48906.1,-56.9652], Tmin=(980.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([O])C[C]=C(18291)',
    structure = SMILES('[CH]=CC([O])C[C]=C'),
    E0 = (606.192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,311.578,311.58,311.585],'cm^-1')),
        HinderedRotor(inertia=(0.00173665,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222757,'amu*angstrom^2'), symmetry=1, barrier=(15.3452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222747,'amu*angstrom^2'), symmetry=1, barrier=(15.3453,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.845989,0.0627885,-5.23917e-05,2.22753e-08,-3.79548e-12,73027,30.3196], Tmin=(100,'K'), Tmax=(1399.61,'K')), NASAPolynomial(coeffs=[14.9752,0.0224084,-9.1157e-06,1.66214e-09,-1.13563e-13,69071.9,-42.596], Tmin=(1399.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=CC([O])[CH]C=C(20301)',
    structure = SMILES('[CH]=CC([O])[CH]C=C'),
    E0 = (438.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,498.898,502.765,503.243],'cm^-1')),
        HinderedRotor(inertia=(0.155342,'amu*angstrom^2'), symmetry=1, barrier=(27.7145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15595,'amu*angstrom^2'), symmetry=1, barrier=(27.7342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156838,'amu*angstrom^2'), symmetry=1, barrier=(27.6954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16974,0.0437017,2.28669e-05,-6.69457e-08,2.99885e-11,52857.2,29.4764], Tmin=(100,'K'), Tmax=(967.874,'K')), NASAPolynomial(coeffs=[18.6372,0.0168291,-5.71725e-06,1.11758e-09,-8.71407e-14,47353.4,-65.1888], Tmin=(967.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH]=[C][CH]C(O)C=C(20060)',
    structure = SMILES('[CH]=[C][CH]C(O)C=C'),
    E0 = (445.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,388.895,389.455],'cm^-1')),
        HinderedRotor(inertia=(0.245568,'amu*angstrom^2'), symmetry=1, barrier=(26.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111417,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244999,'amu*angstrom^2'), symmetry=1, barrier=(26.2217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244898,'amu*angstrom^2'), symmetry=1, barrier=(26.2312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972877,0.0523075,-7.87135e-06,-3.2579e-08,1.73874e-11,53760.1,30.8486], Tmin=(100,'K'), Tmax=(980.407,'K')), NASAPolynomial(coeffs=[17.547,0.0180641,-6.5473e-06,1.24618e-09,-9.27981e-14,48906.1,-56.9652], Tmin=(980.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[CH][O](19109)',
    structure = SMILES('[CH2][C]=C[CH][O]'),
    E0 = (549.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,940.096,1441.46],'cm^-1')),
        HinderedRotor(inertia=(1.73789,'amu*angstrom^2'), symmetry=1, barrier=(39.9575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270964,'amu*angstrom^2'), symmetry=1, barrier=(39.9546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50002,0.0283298,-1.09145e-05,-2.05532e-09,1.67005e-12,66123.7,20.3322], Tmin=(100,'K'), Tmax=(1235.14,'K')), NASAPolynomial(coeffs=[8.27632,0.0176346,-7.65498e-06,1.43667e-09,-9.96445e-14,64085.7,-11.2288], Tmin=(1235.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CCJO) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]=CC([O])=C[CH2](20302)',
    structure = SMILES('[CH2][C]=CC([O])=C[CH2]'),
    E0 = (419.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,373.469,1154.82],'cm^-1')),
        HinderedRotor(inertia=(0.808523,'amu*angstrom^2'), symmetry=1, barrier=(80.325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.808381,'amu*angstrom^2'), symmetry=1, barrier=(80.3227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.808342,'amu*angstrom^2'), symmetry=1, barrier=(80.3102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19426,0.0541726,-2.96037e-05,-4.77445e-09,7.5675e-12,50590.4,26.0115], Tmin=(100,'K'), Tmax=(926.512,'K')), NASAPolynomial(coeffs=[13.4485,0.0205676,-6.44409e-06,1.04412e-09,-6.90838e-14,47491.3,-36.6433], Tmin=(926.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C][CH]C([O])[C]=C(20303)',
    structure = SMILES('C=[C][CH]C([O])[C]=C'),
    E0 = (667.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,562.69,562.845,564.729,564.788],'cm^-1')),
        HinderedRotor(inertia=(0.112337,'amu*angstrom^2'), symmetry=1, barrier=(25.4303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113075,'amu*angstrom^2'), symmetry=1, barrier=(25.4158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113175,'amu*angstrom^2'), symmetry=1, barrier=(25.4255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17167,0.0494391,-1.02245e-05,-2.56226e-08,1.3928e-11,80344.4,30.0143], Tmin=(100,'K'), Tmax=(998.562,'K')), NASAPolynomial(coeffs=[16.1254,0.0184837,-7.20539e-06,1.39075e-09,-1.02723e-13,75914.9,-49.3332], Tmin=(998.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(667.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([O])[CH][C]=C(20304)',
    structure = SMILES('[CH]=CC([O])[CH][C]=C'),
    E0 = (676.335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,500.271,501.05,502.684],'cm^-1')),
        HinderedRotor(inertia=(0.142914,'amu*angstrom^2'), symmetry=1, barrier=(25.5987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145548,'amu*angstrom^2'), symmetry=1, barrier=(25.6375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675589,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13532,0.0483169,-1.97139e-06,-3.76126e-08,1.89947e-11,81460.6,30.0738], Tmin=(100,'K'), Tmax=(980.333,'K')), NASAPolynomial(coeffs=[17.482,0.0162719,-5.96225e-06,1.15894e-09,-8.78379e-14,76590.4,-56.9575], Tmin=(980.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C([O])C=C(20067)',
    structure = SMILES('[CH]=[C][CH]C([O])C=C'),
    E0 = (676.335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,500.652,500.717,500.72],'cm^-1')),
        HinderedRotor(inertia=(0.143985,'amu*angstrom^2'), symmetry=1, barrier=(25.6178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144001,'amu*angstrom^2'), symmetry=1, barrier=(25.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.672559,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13535,0.0483165,-1.96981e-06,-3.76147e-08,1.89956e-11,81460.6,30.0737], Tmin=(100,'K'), Tmax=(980.327,'K')), NASAPolynomial(coeffs=[17.4818,0.0162722,-5.96239e-06,1.15897e-09,-8.78406e-14,76590.5,-56.9568], Tmin=(980.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=CCJC(O)C=C) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C1[CH]CO1(20305)',
    structure = SMILES('C=[C][CH]C1[CH]CO1'),
    E0 = (502.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43065,0.0431282,1.76017e-05,-6.13577e-08,3.00995e-11,60498.6,23.7292], Tmin=(100,'K'), Tmax=(892.965,'K')), NASAPolynomial(coeffs=[14.8997,0.0193273,-3.78594e-06,4.25947e-10,-2.47737e-14,56636.5,-47.8821], Tmin=(892.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(Cds_S) + radical(C=CCJCO) + radical(CCJCO)"""),
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
    label = 'C=C1[CH]C([O])[CH]C1(20229)',
    structure = SMILES('C=C1[CH]C([O])[CH]C1'),
    E0 = (367.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74566,0.0265442,7.20772e-05,-1.15106e-07,4.61696e-11,44343.5,21.6295], Tmin=(100,'K'), Tmax=(964.226,'K')), NASAPolynomial(coeffs=[17.3906,0.0193498,-6.50309e-06,1.29283e-09,-1.02633e-13,38643.8,-67.1893], Tmin=(964.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclopentane) + radical(CCJCO) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=C(O)C=C=C(20306)',
    structure = SMILES('[CH2]C=C(O)C=C=C'),
    E0 = (102.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625392,0.0620265,-2.93286e-05,-1.71872e-08,1.44525e-11,12477.1,23.5744], Tmin=(100,'K'), Tmax=(929.219,'K')), NASAPolynomial(coeffs=[18.7889,0.0150132,-3.76158e-06,5.74945e-10,-4.0105e-14,7755.66,-69.9628], Tmin=(929.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]CC(=O)C=C(18278)',
    structure = SMILES('C=[C]CC(=O)C=C'),
    E0 = (196.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,180,563.298,794.88],'cm^-1')),
        HinderedRotor(inertia=(0.201072,'amu*angstrom^2'), symmetry=1, barrier=(4.62303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0366613,'amu*angstrom^2'), symmetry=1, barrier=(16.4467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.715332,'amu*angstrom^2'), symmetry=1, barrier=(16.4469,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81099,0.0492594,-2.61399e-05,4.85728e-09,-1.65048e-13,23665,23.3833], Tmin=(100,'K'), Tmax=(1906.04,'K')), NASAPolynomial(coeffs=[20.3234,0.0199128,-1.05239e-05,2.01122e-09,-1.34855e-13,14881.6,-82.3979], Tmin=(1906.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]C=C([O])C[CH2](20307)',
    structure = SMILES('[CH2][C]C=C([O])C[CH2]'),
    E0 = (706.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,819.426,4000],'cm^-1')),
        HinderedRotor(inertia=(1.61976,'amu*angstrom^2'), symmetry=1, barrier=(37.2414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0342097,'amu*angstrom^2'), symmetry=1, barrier=(16.3058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00942857,'amu*angstrom^2'), symmetry=1, barrier=(4.49321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709219,'amu*angstrom^2'), symmetry=1, barrier=(16.3063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703097,0.067457,-6.41667e-05,3.19479e-08,-6.34337e-12,85126.7,32.3605], Tmin=(100,'K'), Tmax=(1219.59,'K')), NASAPolynomial(coeffs=[14.5126,0.0221638,-8.45857e-06,1.49551e-09,-1.009e-13,81758.4,-37.0035], Tmin=(1219.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(RCCJ) + radical(C=C(C)OJ) + radical(CCJ2_triplet) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([C]=C)C=O(20137)',
    structure = SMILES('[CH2][CH]C([C]=C)C=O'),
    E0 = (489.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1685,370,1380,1390,370,380,2900,435,360.019,1331.6],'cm^-1')),
        HinderedRotor(inertia=(0.0013008,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0640575,'amu*angstrom^2'), symmetry=1, barrier=(5.87802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0639052,'amu*angstrom^2'), symmetry=1, barrier=(5.87917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.064015,'amu*angstrom^2'), symmetry=1, barrier=(5.88153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44565,0.0606278,-6.17303e-05,3.97758e-08,-1.13266e-11,58926.8,29.7337], Tmin=(100,'K'), Tmax=(824.258,'K')), NASAPolynomial(coeffs=[6.57187,0.0357519,-1.64622e-05,3.16381e-09,-2.22485e-13,58081.7,5.9933], Tmin=(824.258,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S) + radical(CCJCC=O)"""),
)

species(
    label = 'C=[C]C1OC1C=C(20191)',
    structure = SMILES('C=[C]C1OC1C=C'),
    E0 = (305.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41965,0.0414016,2.37626e-05,-7.0931e-08,3.43608e-11,36828.8,24.9581], Tmin=(100,'K'), Tmax=(894.875,'K')), NASAPolynomial(coeffs=[16.6629,0.0149233,-1.68055e-06,4.32877e-11,2.4527e-16,32432.7,-56.2081], Tmin=(894.875,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC=C[CH2](19234)',
    structure = SMILES('[CH2][C]=CC=C[CH2]'),
    E0 = (495.823,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,489.734],'cm^-1')),
        HinderedRotor(inertia=(0.366612,'amu*angstrom^2'), symmetry=1, barrier=(62.3881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366434,'amu*angstrom^2'), symmetry=1, barrier=(62.3853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366714,'amu*angstrom^2'), symmetry=1, barrier=(62.3864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88441,0.0352543,1.66402e-05,-4.67245e-08,2.07763e-11,59720.2,22.7899], Tmin=(100,'K'), Tmax=(947.878,'K')), NASAPolynomial(coeffs=[11.5408,0.0224615,-7.35645e-06,1.26883e-09,-8.84709e-14,56633.7,-29.9045], Tmin=(947.878,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    label = '[CH2][C]=[C]C([O])C=C(20308)',
    structure = SMILES('[CH2][C]=[C]C([O])C=C'),
    E0 = (732.434,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1670,1700,300,440,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0212992,'amu*angstrom^2'), symmetry=1, barrier=(17.0422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231167,'amu*angstrom^2'), symmetry=1, barrier=(17.0414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0212529,'amu*angstrom^2'), symmetry=1, barrier=(17.0323,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30764,0.0595434,-5.2786e-05,2.47038e-08,-4.77445e-12,88188.1,28.408], Tmin=(100,'K'), Tmax=(1214.76,'K')), NASAPolynomial(coeffs=[11.3542,0.0264618,-1.19364e-05,2.28524e-09,-1.60657e-13,85747.3,-22.0152], Tmin=(1214.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(Allyl_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1C[C]=CC1[O](19964)',
    structure = SMILES('[CH2]C1C[C]=CC1[O]'),
    E0 = (509.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73946,0.0327618,4.04088e-05,-7.81237e-08,3.29284e-11,61353.6,24.6873], Tmin=(100,'K'), Tmax=(953.148,'K')), NASAPolynomial(coeffs=[14.85,0.0200836,-6.27374e-06,1.1346e-09,-8.45608e-14,56931,-48.0239], Tmin=(953.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(cyclopentene-vinyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C#CC([O])C=C(20309)',
    structure = SMILES('[CH2]C#CC([O])C=C'),
    E0 = (421.453,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2100,2250,500,550,3010,987.5,1337.5,450,1655,399.422,399.426,399.433],'cm^-1')),
        HinderedRotor(inertia=(0.298113,'amu*angstrom^2'), symmetry=1, barrier=(33.7463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298009,'amu*angstrom^2'), symmetry=1, barrier=(33.7407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00105664,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2711,0.0520796,-2.99208e-05,9.15256e-10,3.61289e-12,50794.3,27.598], Tmin=(100,'K'), Tmax=(1024.49,'K')), NASAPolynomial(coeffs=[12.9248,0.0219512,-8.31542e-06,1.502e-09,-1.04268e-13,47599.8,-32.8434], Tmin=(1024.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CC(C)OJ) + radical(Propargyl)"""),
)

species(
    label = '[CH2]C=[C]C([O])C=C(20310)',
    structure = SMILES('[CH2]C=[C]C([O])C=C'),
    E0 = (494.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,1084.24,1084.3],'cm^-1')),
        HinderedRotor(inertia=(0.818135,'amu*angstrom^2'), symmetry=1, barrier=(18.8105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818154,'amu*angstrom^2'), symmetry=1, barrier=(18.811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818097,'amu*angstrom^2'), symmetry=1, barrier=(18.8097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848085,0.0605163,-4.63971e-05,1.77864e-08,-2.72445e-12,59606.5,29.5987], Tmin=(100,'K'), Tmax=(1547.35,'K')), NASAPolynomial(coeffs=[15.8534,0.0217265,-8.79433e-06,1.58547e-09,-1.06926e-13,54962.8,-49.3439], Tmin=(1547.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=[C]C(O)C=C(20311)',
    structure = SMILES('[CH2][C]=[C]C(O)C=C'),
    E0 = (502.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.171096,'amu*angstrom^2'), symmetry=1, barrier=(3.93384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0399424,'amu*angstrom^2'), symmetry=1, barrier=(15.6744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.039899,'amu*angstrom^2'), symmetry=1, barrier=(15.681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681608,'amu*angstrom^2'), symmetry=1, barrier=(15.6715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1666,0.0632582,-5.76434e-05,2.83305e-08,-5.77743e-12,60486.7,29.108], Tmin=(100,'K'), Tmax=(1154.57,'K')), NASAPolynomial(coeffs=[11.1295,0.0287421,-1.28014e-05,2.43835e-09,-1.71066e-13,58186.1,-20.3895], Tmin=(1154.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=CC([O])[C]=[C]C(20312)',
    structure = SMILES('C=CC([O])[C]=[C]C'),
    E0 = (580.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,249.802,249.83,255.07],'cm^-1')),
        HinderedRotor(inertia=(0.00278138,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255041,'amu*angstrom^2'), symmetry=1, barrier=(10.95,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247292,'amu*angstrom^2'), symmetry=1, barrier=(10.9537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46378,0.0600631,-5.39107e-05,2.7904e-08,-6.301e-12,69958.1,27.578], Tmin=(100,'K'), Tmax=(1019.62,'K')), NASAPolynomial(coeffs=[8.0363,0.0342787,-1.59782e-05,3.10217e-09,-2.1981e-13,68617.8,-4.25831], Tmin=(1019.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C([O])C=[C]C(20313)',
    structure = SMILES('[CH2]C=C([O])C=[C]C'),
    E0 = (301.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771886,0.0647509,-5.86936e-05,2.86869e-08,-5.6272e-12,36405.1,25.9693], Tmin=(100,'K'), Tmax=(1233.45,'K')), NASAPolynomial(coeffs=[13.5925,0.0231748,-8.13312e-06,1.35962e-09,-8.84401e-14,33242.4,-38.5725], Tmin=(1233.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C([O])C=[C]C(20314)',
    structure = SMILES('C=[C]C([O])C=[C]C'),
    E0 = (580.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,249.802,249.83,255.07],'cm^-1')),
        HinderedRotor(inertia=(0.00278138,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255041,'amu*angstrom^2'), symmetry=1, barrier=(10.95,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247292,'amu*angstrom^2'), symmetry=1, barrier=(10.9537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46378,0.0600631,-5.39107e-05,2.7904e-08,-6.301e-12,69958.1,27.578], Tmin=(100,'K'), Tmax=(1019.62,'K')), NASAPolynomial(coeffs=[8.0363,0.0342787,-1.59782e-05,3.10217e-09,-2.1981e-13,68617.8,-4.25831], Tmin=(1019.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([O])C=[C]C(20315)',
    structure = SMILES('[CH]=CC([O])C=[C]C'),
    E0 = (590.19,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,320.041,320.042],'cm^-1')),
        HinderedRotor(inertia=(0.155467,'amu*angstrom^2'), symmetry=1, barrier=(11.2997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155462,'amu*angstrom^2'), symmetry=1, barrier=(11.2998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155462,'amu*angstrom^2'), symmetry=1, barrier=(11.2999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3323,0.0602,-5.05887e-05,2.28415e-08,-4.33932e-12,71078.1,27.9682], Tmin=(100,'K'), Tmax=(1216.54,'K')), NASAPolynomial(coeffs=[10.433,0.0302766,-1.36928e-05,2.62235e-09,-1.84237e-13,68863.8,-17.721], Tmin=(1216.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[O]C1[CH]CC[C]=C1(20038)',
    structure = SMILES('[O]C1[CH]CC[C]=C1'),
    E0 = (474.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01048,0.0236821,6.67032e-05,-1.03004e-07,4.05848e-11,57102.2,23.7027], Tmin=(100,'K'), Tmax=(969.931,'K')), NASAPolynomial(coeffs=[14.725,0.021493,-7.61654e-06,1.48813e-09,-1.14281e-13,52272.3,-49.4332], Tmin=(969.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CCJCO) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=C([O])C=C(20316)',
    structure = SMILES('C=CC=C([O])C=C'),
    E0 = (69.6051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689397,0.0604166,-2.702e-05,-1.94489e-08,1.53256e-11,8502.2,23.2573], Tmin=(100,'K'), Tmax=(927.426,'K')), NASAPolynomial(coeffs=[18.89,0.0137934,-3.1686e-06,4.66226e-10,-3.28656e-14,3755.37,-70.5697], Tmin=(927.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.6051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C]=[C]C([O])C[CH2](20317)',
    structure = SMILES('[CH2][C]=[C]C([O])C[CH2]'),
    E0 = (813.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,290.885,290.892,1390.19],'cm^-1')),
        HinderedRotor(inertia=(0.2955,'amu*angstrom^2'), symmetry=1, barrier=(17.7431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166437,'amu*angstrom^2'), symmetry=1, barrier=(9.99291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0129376,'amu*angstrom^2'), symmetry=1, barrier=(17.7432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0129376,'amu*angstrom^2'), symmetry=1, barrier=(17.7431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.046,0.0656527,-6.13578e-05,3.08423e-08,-6.3968e-12,97965.9,30.8244], Tmin=(100,'K'), Tmax=(1141.37,'K')), NASAPolynomial(coeffs=[11.6157,0.0286107,-1.26769e-05,2.40816e-09,-1.68745e-13,95553.2,-21.5659], Tmin=(1141.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(813.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([O])[CH]C(20318)',
    structure = SMILES('[CH2][C]=[C]C([O])[CH]C'),
    E0 = (808.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,1670,1700,300,440,3025,407.5,1350,352.5,180,180,784.066],'cm^-1')),
        HinderedRotor(inertia=(0.148258,'amu*angstrom^2'), symmetry=1, barrier=(3.40875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00781376,'amu*angstrom^2'), symmetry=1, barrier=(3.40882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0384176,'amu*angstrom^2'), symmetry=1, barrier=(16.7595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0384172,'amu*angstrom^2'), symmetry=1, barrier=(16.7596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928441,0.0629491,-5.3574e-05,2.35797e-08,-4.19334e-12,97331.7,31.9549], Tmin=(100,'K'), Tmax=(1336.61,'K')), NASAPolynomial(coeffs=[13.9105,0.0240987,-9.97464e-06,1.83357e-09,-1.25961e-13,93861.3,-34.4425], Tmin=(1336.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P) + radical(Cds_S) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC1[CH]C(=C)O1(20155)',
    structure = SMILES('C=CC1[CH]C(=C)O1'),
    E0 = (82.1486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62324,0.0176952,0.000124065,-1.94081e-07,8.18468e-11,9997.98,20.4662], Tmin=(100,'K'), Tmax=(916.225,'K')), NASAPolynomial(coeffs=[26.1871,0.000159011,5.91666e-06,-1.25618e-09,7.61496e-14,1731.63,-116.438], Tmin=(916.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.1486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=CCJC(O)C=C)"""),
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
    E0 = (429.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (537.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (582.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (524.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (473.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (663.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (552.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (559.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (527.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (543.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (587.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (593.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (578.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (738.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (554.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (751.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (643.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (490.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (650.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (747.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (479.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (838.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (631.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (878.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (888.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (888.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (561.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (581.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (472.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (452.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (452.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (729.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (652.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (432.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (902.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (739.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (944.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (549.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (648.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (441.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (698.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (625.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (742.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (560.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (723.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (722.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (516.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (492.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (877.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (833.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (437.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['C=CC=O(5269)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2]C1OC1[CH][C]=C(20294)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2]C1C([O])C1[C]=C(20134)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(153.076,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 152.9 to 153.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2]C1C(=C)[CH]C1[O](20100)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.16717e+10,'s^-1'), n=0.521143, Ea=(95.4406,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra_2H_pri;radadd_intra] for rate rule [R5;doublebond_intra_2H_pri;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH2]C=C([O])C=C=C(20295)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', '[CH]=C=CC([O])C=C(20057)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=CC=O(5269)', '[CH][C]=C(18825)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.5e+06,'cm^3/(mol*s)'), n=2.16, Ea=(19.2464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdH_O;YJ] for rate rule [CO-CdH_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(T)(63)', '[CH2]C=CC=C=C(19683)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(64)', 'C=[C]C=C[O](18052)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0131003,'m^3/(mol*s)'), n=2.40999, Ea=(12.7705,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CdsJ-H] for rate rule [CO_O;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2][C]=CC(O)=C[CH2](20296)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2]C=C([O])C[C]=C(18288)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C[CH]C([O])C=C(20059)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.3587e+07,'s^-1'), n=1.74167, Ea=(154.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]CC([O])[C]=C(18290)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2]C=CC([O])=C[CH2](20298)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.30234e+06,'s^-1'), n=1.68744, Ea=(125.264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]CC([O])C=C(18292)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C([O])[CH]C=C(20299)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_doubleC] for rate rule [R4HJ_2;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC(O)[CH][C]=C(20300)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CC([O])C[C]=C(18291)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC([O])[CH]C=C(20301)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C][CH]C(O)C=C(20060)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(64)', '[CH2][C]=C[CH][O](19109)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.88428e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2][C]=CC([O])=C[CH2](20302)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', 'C=[C][CH]C([O])[C]=C(20303)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH]=CC([O])[CH][C]=C(20304)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH]=[C][CH]C([O])C=C(20067)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['C=[C][CH]C1[CH]CO1(20305)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.15968e+08,'s^-1'), n=1.10215, Ea=(132.51,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['C=[C]C1C[CH]C1[O](20262)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.61918e+08,'s^-1'), n=0.930343, Ea=(152.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 150.6 to 152.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['C=C1[CH]C([O])[CH]C1(20229)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.82632e+08,'s^-1'), n=0.716884, Ea=(43.4606,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2]C=C(O)C=C=C(20306)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['C=[C]CC(=O)C=C(18278)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]C=C([O])C[CH2](20307)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C([C]=C)C=O(20137)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['C=[C]C1OC1C=C(20191)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(T)(63)', '[CH2][C]=CC=C[CH2](19234)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C=C[O](5266)', '[CH][C]=C(18825)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(8)', '[CH2][C]=[C]C([O])C=C(20308)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2]C1C[C]=CC1[O](19964)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(9.40382e+09,'s^-1'), n=0.352, Ea=(120.148,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(8)', '[CH2]C#CC([O])C=C(20309)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C=C[O](5266)', 'C#C[CH2](17441)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C=[C]C([O])C=C(20310)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][C]=[C]C(O)C=C(20311)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(117344,'s^-1'), n=2.01217, Ea=(123.77,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;O_H_out] + [R3H_SS_Cs;Cd_rad_out;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=CC([O])[C]=[C]C(20312)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[CH2]C=C([O])C=[C]C(20313)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=[C]C([O])C=[C]C(20314)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=CC([O])C=[C]C(20315)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_4;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['[O]C1[CH]CC[C]=C1(20038)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['C=CC=C([O])C=C(20316)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][C]=[C]C([O])C[CH2](20317)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2][C]=[C]C([O])[CH]C(20318)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=[C][CH]C([O])C=C(18289)'],
    products = ['C=CC1[CH]C(=C)O1(20155)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

network(
    label = '4188',
    isomers = [
        'C=[C][CH]C([O])C=C(18289)',
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
    label = '4188',
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

