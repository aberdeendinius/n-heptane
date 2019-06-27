species(
    label = '[CH2]C([CH2])C=C=C[O](22494)',
    structure = SMILES('[CH2]C([CH2])C=C=C[O]'),
    E0 = (435.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,366.807,366.808],'cm^-1')),
        HinderedRotor(inertia=(0.147672,'amu*angstrom^2'), symmetry=1, barrier=(14.0995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20092,'amu*angstrom^2'), symmetry=1, barrier=(114.663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20092,'amu*angstrom^2'), symmetry=1, barrier=(114.662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.164353,0.071034,-6.61192e-05,3.11009e-08,-5.50703e-12,52502.8,32.0379], Tmin=(100,'K'), Tmax=(1616.63,'K')), NASAPolynomial(coeffs=[17.7979,0.0144447,-2.34316e-06,1.53703e-10,-2.60997e-15,48282.2,-58.339], Tmin=(1616.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(Isobutyl) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=C(87)',
    structure = SMILES('[CH2]C=C'),
    E0 = (157.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.570287,'amu*angstrom^2'), symmetry=1, barrier=(32.8573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3193,0.00566487,4.27449e-05,-5.78831e-08,2.21699e-11,18990.6,9.19646], Tmin=(100,'K'), Tmax=(951.999,'K')), NASAPolynomial(coeffs=[7.55715,0.0114811,-3.63952e-06,6.63584e-10,-4.95318e-14,17113.3,-16.6624], Tmin=(951.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
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
    label = '[CH2]C1CC1[C]=C[O](23586)',
    structure = SMILES('[CH2]C1CC1[C]=C[O]'),
    E0 = (457.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30956,0.0391041,3.83004e-05,-8.96923e-08,4.10211e-11,55091.3,26.0847], Tmin=(100,'K'), Tmax=(920.655,'K')), NASAPolynomial(coeffs=[19.7318,0.0110888,-8.17275e-07,1.19716e-11,-4.75583e-15,49494.4,-73.2433], Tmin=(920.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(Cds_S) + radical(Isobutyl) + radical(C=COJ)"""),
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
    label = '[CH2]C(=C)C=C=C[O](22603)',
    structure = SMILES('[CH2]C(=C)C=C=C[O]'),
    E0 = (281.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43154,'amu*angstrom^2'), symmetry=1, barrier=(32.914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4384,'amu*angstrom^2'), symmetry=1, barrier=(33.0715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612281,0.0583592,-1.49787e-05,-3.75866e-08,2.29716e-11,33942.8,22.4721], Tmin=(100,'K'), Tmax=(929.873,'K')), NASAPolynomial(coeffs=[21.882,0.0082147,-7.92857e-07,6.55313e-11,-8.55619e-15,28199.5,-88.2081], Tmin=(929.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C=C=C=O(23587)',
    structure = SMILES('[CH2]C([CH2])C=C=C=O'),
    E0 = (479.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,540,610,2055,1380,1390,370,380,2900,435,2120,512.5,787.5,3010,987.5,1337.5,450,1655,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000462027,'amu*angstrom^2'), symmetry=1, barrier=(5.24588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000458271,'amu*angstrom^2'), symmetry=1, barrier=(5.20324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000459508,'amu*angstrom^2'), symmetry=1, barrier=(5.21728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0741,0.0363494,-9.63264e-06,-1.55869e-08,9.9347e-12,57713.7,11.1927], Tmin=(100,'K'), Tmax=(906.569,'K')), NASAPolynomial(coeffs=[10.1392,0.0170984,-4.80616e-06,7.38044e-10,-4.77792e-14,55580.2,-30.6275], Tmin=(906.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = 'C=CC=C=C[O](22346)',
    structure = SMILES('C=CC=C=C[O]'),
    E0 = (167.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40235,'amu*angstrom^2'), symmetry=1, barrier=(32.2429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43985,0.0410476,6.30691e-06,-5.03167e-08,2.59189e-11,20212.2,19.3136], Tmin=(100,'K'), Tmax=(926.458,'K')), NASAPolynomial(coeffs=[18.7967,0.00545165,2.41126e-07,-1.15506e-10,3.7087e-15,15307.7,-72.2094], Tmin=(926.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
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
    label = '[CH2]C(C)=C[C]=C[O](22744)',
    structure = SMILES('[CH2]C(C)=C[C]=C[O]'),
    E0 = (268.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.68877,'amu*angstrom^2'), symmetry=1, barrier=(38.8282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69537,'amu*angstrom^2'), symmetry=1, barrier=(38.9799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.70347,'amu*angstrom^2'), symmetry=1, barrier=(39.1662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769871,0.0589104,-2.05177e-05,-2.54944e-08,1.74411e-11,32428.1,24.5766], Tmin=(100,'K'), Tmax=(916.919,'K')), NASAPolynomial(coeffs=[17.731,0.016892,-4.08485e-06,5.87727e-10,-3.91937e-14,27973.6,-63.109], Tmin=(916.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C([CH2])[CH]C#CO(23588)',
    structure = SMILES('[CH2]C([CH2])[CH]C#CO'),
    E0 = (507.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,2100,2250,500,550,1380,1390,370,380,2900,435,1768.31,3090.42],'cm^-1')),
        HinderedRotor(inertia=(1.82075,'amu*angstrom^2'), symmetry=1, barrier=(77.0181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82075,'amu*angstrom^2'), symmetry=1, barrier=(77.0181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.011364,'amu*angstrom^2'), symmetry=1, barrier=(77.0181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82076,'amu*angstrom^2'), symmetry=1, barrier=(77.0181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82075,'amu*angstrom^2'), symmetry=1, barrier=(77.0181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771532,0.0671463,-7.16884e-05,4.33221e-08,-1.02975e-11,61204.5,29.7572], Tmin=(100,'K'), Tmax=(1143.79,'K')), NASAPolynomial(coeffs=[11.6571,0.0232571,-6.49693e-06,8.75254e-10,-4.73557e-14,59095.1,-22.5571], Tmin=(1143.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Sec_Propargyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)[C]=C=C[O](22745)',
    structure = SMILES('[CH2]C(C)[C]=C=C[O]'),
    E0 = (467.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,269.233,270.589],'cm^-1')),
        HinderedRotor(inertia=(0.372214,'amu*angstrom^2'), symmetry=1, barrier=(19.5179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.380053,'amu*angstrom^2'), symmetry=1, barrier=(19.5186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374958,'amu*angstrom^2'), symmetry=1, barrier=(19.5142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.727262,0.0630504,-4.11983e-05,7.30549e-10,6.8282e-12,56401.3,28.2086], Tmin=(100,'K'), Tmax=(934.535,'K')), NASAPolynomial(coeffs=[16.3486,0.0185262,-5.58841e-06,9.05341e-10,-6.09109e-14,52506,-51.3169], Tmin=(934.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH2])C#C[CH]O(23589)',
    structure = SMILES('[CH2]C([CH2])C#C[CH]O'),
    E0 = (511.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,2100,2250,500,550,1380,1390,370,380,2900,435,180,1385.51],'cm^-1')),
        HinderedRotor(inertia=(0.145168,'amu*angstrom^2'), symmetry=1, barrier=(3.3377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.78659,'amu*angstrom^2'), symmetry=1, barrier=(64.0691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145129,'amu*angstrom^2'), symmetry=1, barrier=(3.3368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00244821,'amu*angstrom^2'), symmetry=1, barrier=(3.33496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.78654,'amu*angstrom^2'), symmetry=1, barrier=(64.0679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16993,0.0653574,-5.73354e-05,5.51988e-09,1.64553e-11,61657.2,29.6981], Tmin=(100,'K'), Tmax=(654.169,'K')), NASAPolynomial(coeffs=[11.6143,0.0220715,-5.26566e-06,5.41154e-10,-1.8678e-14,59850.4,-19.6229], Tmin=(654.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCsJOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)[CH][C]=C=O(22746)',
    structure = SMILES('[CH2]C(C)[CH][C]=C=O'),
    E0 = (396.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1685,370,2120,512.5,787.5,271.73,3179.62],'cm^-1')),
        HinderedRotor(inertia=(0.00228291,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.688664,'amu*angstrom^2'), symmetry=1, barrier=(36.0851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45949,'amu*angstrom^2'), symmetry=1, barrier=(76.4757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45946,'amu*angstrom^2'), symmetry=1, barrier=(76.4751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13353,0.0640608,-6.05561e-05,3.28352e-08,-7.47043e-12,47775.9,27.5445], Tmin=(100,'K'), Tmax=(1042.07,'K')), NASAPolynomial(coeffs=[9.75022,0.0309853,-1.29456e-05,2.37614e-09,-1.63051e-13,45980,-14.3811], Tmin=(1042.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCJC(C)=C=O) + radical(Isobutyl) + radical(CCCJ=C=O)"""),
)

species(
    label = '[CH2]C([CH2])=C[C]=CO(23590)',
    structure = SMILES('[CH2]C([CH2])=C[C]=CO'),
    E0 = (245.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.641515,0.0568858,-8.4372e-09,-5.91579e-08,3.32038e-11,29622.3,25.1882], Tmin=(100,'K'), Tmax=(890.96,'K')), NASAPolynomial(coeffs=[21.8566,0.00890798,1.18661e-06,-5.06436e-10,3.8169e-14,23965.9,-85.2412], Tmin=(890.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][CH][CH2](6136)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00132719,'amu*angstrom^2'), symmetry=1, barrier=(2.41051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00132749,'amu*angstrom^2'), symmetry=1, barrier=(2.41088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34321,0.013802,2.16426e-06,-5.76329e-09,1.61332e-12,58271,14.955], Tmin=(100,'K'), Tmax=(1447.11,'K')), NASAPolynomial(coeffs=[4.39505,0.0167645,-6.99091e-06,1.25741e-09,-8.38108e-14,57351.9,7.36811], Tmin=(1447.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])=C[C]=C[O](23591)',
    structure = SMILES('[CH2]C([CH2])=C[C]=C[O]'),
    E0 = (386.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,211.904,211.962],'cm^-1')),
        HinderedRotor(inertia=(1.55124,'amu*angstrom^2'), symmetry=1, barrier=(49.5142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5558,'amu*angstrom^2'), symmetry=1, barrier=(49.5206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55615,'amu*angstrom^2'), symmetry=1, barrier=(49.5177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05505,0.0499202,3.21203e-06,-5.23893e-08,2.80469e-11,46619.3,25.1123], Tmin=(100,'K'), Tmax=(903.846,'K')), NASAPolynomial(coeffs=[18.4767,0.0128071,-1.55784e-06,7.6688e-11,-3.76707e-15,41836.6,-66.2116], Tmin=(903.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([CH2])[C]=C=C[O](23592)',
    structure = SMILES('[CH2]C([CH2])[C]=C=C[O]'),
    E0 = (672.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,540,610,2055,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,368.949,372.833],'cm^-1')),
        HinderedRotor(inertia=(0.113988,'amu*angstrom^2'), symmetry=1, barrier=(11.1672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22919,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19802,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280306,0.0702613,-7.32263e-05,3.88442e-08,-7.83627e-12,81084.8,30.8976], Tmin=(100,'K'), Tmax=(1368.07,'K')), NASAPolynomial(coeffs=[16.8728,0.0138306,-2.67298e-06,2.33087e-10,-7.48051e-15,77285.7,-51.6436], Tmin=(1368.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Isobutyl) + radical(Isobutyl) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH2])[CH][C]=C=O(23593)',
    structure = SMILES('[CH2]C([CH2])[CH][C]=C=O'),
    E0 = (601.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2120,512.5,787.5,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,353.873,354.029],'cm^-1')),
        HinderedRotor(inertia=(0.00134601,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134535,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134613,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625164,'amu*angstrom^2'), symmetry=1, barrier=(55.5657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12979,0.0657864,-7.198e-05,4.16895e-08,-8.22884e-12,72440.3,28.6589], Tmin=(100,'K'), Tmax=(730.007,'K')), NASAPolynomial(coeffs=[9.71738,0.027342,-1.06778e-05,1.86362e-09,-1.23276e-13,70957,-11.6402], Tmin=(730.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(Isobutyl) + radical(CCJC(C)=C=O) + radical(Isobutyl) + radical(CCCJ=C=O)"""),
)

species(
    label = '[CH2]C([CH2])C=C1[CH]O1(23594)',
    structure = SMILES('[CH2]C([CH2])C=C1[CH]O1'),
    E0 = (461.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17039,0.0392327,4.7621e-05,-1.08208e-07,5.0349e-11,55607.6,25.8168], Tmin=(100,'K'), Tmax=(904.679,'K')), NASAPolynomial(coeffs=[22.7199,0.00501257,3.11885e-06,-8.08546e-10,5.37332e-14,49209.8,-89.7986], Tmin=(904.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(methyleneoxirane) + radical(C=CCJO) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[CH]C(=C[O])C1(23595)',
    structure = SMILES('[CH2]C1[CH]C(=C[O])C1'),
    E0 = (346.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.99481,-0.0138987,0.000130825,-1.27166e-07,2.97288e-11,41342.8,-15.151], Tmin=(100,'K'), Tmax=(1686.7,'K')), NASAPolynomial(coeffs=[70.5094,0.0272531,-6.99946e-05,1.71171e-08,-1.27729e-12,-5338.83,-416.602], Tmin=(1686.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(methylenecyclobutane) + radical(Allyl_S) + radical(Isobutyl) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH2])C1[C]=CO1(23596)',
    structure = SMILES('[CH2]C([CH2])C1[C]=CO1'),
    E0 = (553.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.950074,0.044873,3.43219e-05,-9.57725e-08,4.61111e-11,66750.3,26.2796], Tmin=(100,'K'), Tmax=(906.052,'K')), NASAPolynomial(coeffs=[23.4051,0.00478318,2.94314e-06,-7.61385e-10,5.03553e-14,60257.7,-93.2118], Tmin=(906.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1C=[C]C([O])C1(23597)',
    structure = SMILES('[CH2]C1C=[C]C([O])C1'),
    E0 = (510.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73772,0.0325011,4.18422e-05,-8.02504e-08,3.38744e-11,61482.6,24.7282], Tmin=(100,'K'), Tmax=(950.687,'K')), NASAPolynomial(coeffs=[15.0722,0.0195976,-5.96173e-06,1.07116e-09,-8.0153e-14,56995,-49.1962], Tmin=(950.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(C)C=C=C[O](22735)',
    structure = SMILES('C=C(C)C=C=C[O]'),
    E0 = (129.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18958,'amu*angstrom^2'), symmetry=1, barrier=(27.3509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18931,'amu*angstrom^2'), symmetry=1, barrier=(27.3445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.602559,0.0609861,-2.40908e-05,-2.33147e-08,1.64695e-11,15719.7,22.919], Tmin=(100,'K'), Tmax=(942.614,'K')), NASAPolynomial(coeffs=[19.6721,0.0141266,-3.7268e-06,6.1906e-10,-4.57526e-14,10611.4,-75.9807], Tmin=(942.614,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(=C)C=C=CO(23598)',
    structure = SMILES('[CH2]C(=C)C=C=CO'),
    E0 = (139.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.198633,0.0653397,-1.83301e-05,-4.40387e-08,2.79111e-11,16945.8,22.5476], Tmin=(100,'K'), Tmax=(913.206,'K')), NASAPolynomial(coeffs=[25.212,0.004399,1.90399e-06,-5.06434e-10,3.24597e-14,10349.9,-106.957], Tmin=(913.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)C=C=C=O(22736)',
    structure = SMILES('[CH2]C(C)C=C=C=O'),
    E0 = (274.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000680058,'amu*angstrom^2'), symmetry=1, barrier=(7.72142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000678588,'amu*angstrom^2'), symmetry=1, barrier=(7.70473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000680926,'amu*angstrom^2'), symmetry=1, barrier=(7.73128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02709,0.035033,1.42948e-06,-2.59453e-08,1.24818e-11,33051.6,10.2713], Tmin=(100,'K'), Tmax=(966.365,'K')), NASAPolynomial(coeffs=[10.4,0.0203801,-6.87709e-06,1.20619e-09,-8.40114e-14,30499.3,-34.6693], Tmin=(966.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C[C]=C[O](23599)',
    structure = SMILES('[CH2][C]([CH2])C[C]=C[O]'),
    E0 = (692.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,360,370,350,180,742.276,743.408],'cm^-1')),
        HinderedRotor(inertia=(0.121155,'amu*angstrom^2'), symmetry=1, barrier=(2.78559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120427,'amu*angstrom^2'), symmetry=1, barrier=(2.76885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121017,'amu*angstrom^2'), symmetry=1, barrier=(2.78242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121934,'amu*angstrom^2'), symmetry=1, barrier=(2.8035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01075,0.0689468,-7.29454e-05,3.62732e-08,-3.3278e-12,83450.6,30.7939], Tmin=(100,'K'), Tmax=(696.827,'K')), NASAPolynomial(coeffs=[10.4388,0.0272674,-1.00052e-05,1.67778e-09,-1.07933e-13,81834.6,-13.4526], Tmin=(696.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(692.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Isobutyl) + radical(Cds_S) + radical(Tertalkyl) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[C]C=C[O](23600)',
    structure = SMILES('[CH2]C([CH2])[C]C=C[O]'),
    E0 = (717.926,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,413.974,414.536,414.855,414.981],'cm^-1')),
        HinderedRotor(inertia=(0.0108666,'amu*angstrom^2'), symmetry=1, barrier=(1.32062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149385,'amu*angstrom^2'), symmetry=1, barrier=(18.2423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623655,'amu*angstrom^2'), symmetry=1, barrier=(76.1302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.624325,'amu*angstrom^2'), symmetry=1, barrier=(76.1247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.230013,0.0737228,-7.10598e-05,3.43757e-08,-6.25732e-12,86515.4,33.5237], Tmin=(100,'K'), Tmax=(1562.55,'K')), NASAPolynomial(coeffs=[18.5715,0.0137526,-2.12422e-06,1.14751e-10,1.39351e-16,82085.2,-60.9489], Tmin=(1562.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(717.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(Isobutyl) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH2])[CH]C=[C][O](23601)',
    structure = SMILES('[CH2]C([CH2])[CH]C=[C][O]'),
    E0 = (650.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,391.551,392.211,392.305],'cm^-1')),
        HinderedRotor(inertia=(0.00109529,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00110057,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00109361,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00109838,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1022,0.0570387,-3.43875e-05,-4.26534e-11,5.95561e-12,78355.7,30.851], Tmin=(100,'K'), Tmax=(918.393,'K')), NASAPolynomial(coeffs=[13.0476,0.0226897,-7.15974e-06,1.15236e-09,-7.52596e-14,75416.1,-29.8207], Tmin=(918.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(650.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Isobutyl) + radical(C=CJO) + radical(Isobutyl) + radical(Allyl_S) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH2])=C[CH][CH][O](23602)',
    structure = SMILES('[CH2]C([CH2])=C[CH][CH][O]'),
    E0 = (574.488,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,1273.88],'cm^-1')),
        HinderedRotor(inertia=(1.10237,'amu*angstrom^2'), symmetry=1, barrier=(25.3457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.405542,'amu*angstrom^2'), symmetry=1, barrier=(55.8946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115824,'amu*angstrom^2'), symmetry=1, barrier=(55.8949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115825,'amu*angstrom^2'), symmetry=1, barrier=(55.8951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.758659,0.0750419,-7.94836e-05,4.67425e-08,-1.14304e-11,69208.6,25.1628], Tmin=(100,'K'), Tmax=(974.976,'K')), NASAPolynomial(coeffs=[10.8976,0.033446,-1.54896e-05,2.98565e-09,-2.10577e-13,67231.5,-23.4949], Tmin=(974.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCsJOH) + radical(Allyl_P) + radical(C=CCJCO) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([CH2])=C[C]C[O](23603)',
    structure = SMILES('[CH2]C([CH2])=C[C]C[O]'),
    E0 = (725.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,229.595,229.596,917.19,3934.94],'cm^-1')),
        HinderedRotor(inertia=(0.530902,'amu*angstrom^2'), symmetry=1, barrier=(19.8592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34578,'amu*angstrom^2'), symmetry=1, barrier=(50.3417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34577,'amu*angstrom^2'), symmetry=1, barrier=(50.3417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0477692,'amu*angstrom^2'), symmetry=1, barrier=(50.3417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1585,0.0658697,-6.19014e-05,3.29838e-08,-7.49089e-12,87355.7,26.8815], Tmin=(100,'K'), Tmax=(1030.12,'K')), NASAPolynomial(coeffs=[9.35726,0.0340337,-1.55439e-05,2.98253e-09,-2.09904e-13,85666.5,-12.9161], Tmin=(1030.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(725.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCJ2_triplet) + radical(Allyl_P) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = '[CH2][CH]CC=C=C[O](22491)',
    structure = SMILES('[CH2][CH]CC=C=C[O]'),
    E0 = (432.639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,815.652,816.721,817.009],'cm^-1')),
        HinderedRotor(inertia=(0.384228,'amu*angstrom^2'), symmetry=1, barrier=(9.87038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428246,'amu*angstrom^2'), symmetry=1, barrier=(9.84621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26317,'amu*angstrom^2'), symmetry=1, barrier=(29.0427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773252,0.0626786,-5.25105e-05,2.29499e-08,-4.01481e-12,52157.6,30.5618], Tmin=(100,'K'), Tmax=(1373.84,'K')), NASAPolynomial(coeffs=[14.7475,0.0219923,-8.08813e-06,1.39378e-09,-9.22196e-14,48317.9,-41.2942], Tmin=(1373.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(RCCJC) + radical(RCCJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CC1CC1(23604)',
    structure = SMILES('[O]C=C=CC1CC1'),
    E0 = (187.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37,0.0355618,5.04473e-05,-1.01084e-07,4.4069e-11,22674.7,23.1996], Tmin=(100,'K'), Tmax=(938.713,'K')), NASAPolynomial(coeffs=[20.208,0.0118377,-2.00214e-06,3.37694e-10,-3.22135e-14,16646.6,-79.7619], Tmin=(938.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=C[C]=C[O](22677)',
    structure = SMILES('[CH2]C=C[C]=C[O]'),
    E0 = (307.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.88234,'amu*angstrom^2'), symmetry=1, barrier=(43.2786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87521,'amu*angstrom^2'), symmetry=1, barrier=(43.1147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60086,0.0390403,9.67809e-06,-5.23095e-08,2.68547e-11,37097,20.8382], Tmin=(100,'K'), Tmax=(904.514,'K')), NASAPolynomial(coeffs=[16.9126,0.00812039,-6.13762e-08,-1.59928e-10,1.13511e-14,32822,-59.8156], Tmin=(904.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
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
    label = 'C#C[CH]C([CH2])[CH2](19176)',
    structure = SMILES('C#C[CH]C([CH2])[CH2]'),
    E0 = (649.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,937.332],'cm^-1')),
        HinderedRotor(inertia=(0.220286,'amu*angstrom^2'), symmetry=1, barrier=(5.06482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.099972,'amu*angstrom^2'), symmetry=1, barrier=(62.1928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221828,'amu*angstrom^2'), symmetry=1, barrier=(5.10025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.099932,'amu*angstrom^2'), symmetry=1, barrier=(62.1682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0061,0.0596319,-6.11782e-05,3.54797e-08,-7.92933e-12,78212.5,25.573], Tmin=(100,'K'), Tmax=(1277.63,'K')), NASAPolynomial(coeffs=[10.8808,0.0199365,-4.26605e-06,4.0436e-10,-1.34998e-14,76405.8,-21.6817], Tmin=(1277.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C([CH2])C=C=C[O](23605)',
    structure = SMILES('[CH]C([CH2])C=C=C[O]'),
    E0 = (678.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,442.757,442.77,442.778,442.811,442.839],'cm^-1')),
        HinderedRotor(inertia=(0.147788,'amu*angstrom^2'), symmetry=1, barrier=(20.5621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000859779,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622377,'amu*angstrom^2'), symmetry=1, barrier=(86.5883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.13345,0.0711944,-6.88209e-05,3.25614e-08,-5.78154e-12,81743.1,31.8436], Tmin=(100,'K'), Tmax=(1592.47,'K')), NASAPolynomial(coeffs=[19.4027,0.0110523,-1.74312e-06,1.14591e-10,-2.42255e-15,76924.7,-67.0895], Tmin=(1592.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(678.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH2])C#CC=O(23606)',
    structure = SMILES('[CH2]C([CH2])C#CC=O'),
    E0 = (396.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2100,2250,500,550,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1973.52],'cm^-1')),
        HinderedRotor(inertia=(0.550846,'amu*angstrom^2'), symmetry=1, barrier=(12.665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.05471,'amu*angstrom^2'), symmetry=1, barrier=(70.2337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0254445,'amu*angstrom^2'), symmetry=1, barrier=(70.2468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.05258,'amu*angstrom^2'), symmetry=1, barrier=(70.1849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38395,0.0565302,-4.88106e-05,1.97067e-08,-1.21168e-12,47828.9,26.9809], Tmin=(100,'K'), Tmax=(817.046,'K')), NASAPolynomial(coeffs=[9.89746,0.0246559,-8.2943e-06,1.33549e-09,-8.4686e-14,46110.5,-14.3744], Tmin=(817.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[C]=CC=O(23607)',
    structure = SMILES('[CH2]C([CH2])[C]=CC=O'),
    E0 = (476.244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,196.675],'cm^-1')),
        HinderedRotor(inertia=(0.372849,'amu*angstrom^2'), symmetry=1, barrier=(10.2246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00435938,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0036532,'amu*angstrom^2'), symmetry=1, barrier=(10.2244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.39762,'amu*angstrom^2'), symmetry=1, barrier=(65.8256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22161,0.0647606,-6.97008e-05,4.58043e-08,-1.28469e-11,57375.8,28.4827], Tmin=(100,'K'), Tmax=(849.793,'K')), NASAPolynomial(coeffs=[7.82074,0.0336988,-1.48735e-05,2.79275e-09,-1.93552e-13,56254.2,-2.28025], Tmin=(849.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Isobutyl) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH]C=C=O(23608)',
    structure = SMILES('[CH2]C([CH2])[CH]C=C=O'),
    E0 = (399.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04476,0.0661429,-7.08581e-05,4.5038e-08,-1.18533e-11,48115.8,29.2415], Tmin=(100,'K'), Tmax=(917.052,'K')), NASAPolynomial(coeffs=[9.33077,0.0300006,-1.17403e-05,2.06071e-09,-1.36984e-13,46596.1,-10.0162], Tmin=(917.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(Isobutyl) + radical(CCJC(C)=C=O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])=CC=C[O](18261)',
    structure = SMILES('[CH2]C([CH2])=CC=C[O]'),
    E0 = (187.622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,302.144,304.219],'cm^-1')),
        HinderedRotor(inertia=(0.705847,'amu*angstrom^2'), symmetry=1, barrier=(46.0574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696885,'amu*angstrom^2'), symmetry=1, barrier=(46.0341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.700917,'amu*angstrom^2'), symmetry=1, barrier=(46.0506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12473,0.0434961,3.17733e-05,-8.43769e-08,3.93298e-11,22687.6,25.2642], Tmin=(100,'K'), Tmax=(920.678,'K')), NASAPolynomial(coeffs=[19.9679,0.0128001,-1.58415e-06,1.4475e-10,-1.32279e-14,17049.2,-75.8635], Tmin=(920.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C1[CH][C]=COC1(23609)',
    structure = SMILES('[CH2]C1[CH][C]=COC1'),
    E0 = (396.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81608,0.0269297,6.88301e-05,-1.18905e-07,5.15188e-11,47825.8,18.4541], Tmin=(100,'K'), Tmax=(904.945,'K')), NASAPolynomial(coeffs=[17.4509,0.0132697,-4.36741e-07,-1.68005e-10,1.12289e-14,42725.7,-67.9579], Tmin=(904.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(Allyl_S) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=C)C=CC=O(18274)',
    structure = SMILES('[CH2]C(=C)C=CC=O'),
    E0 = (84.3418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.962873,'amu*angstrom^2'), symmetry=1, barrier=(22.1383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.963335,'amu*angstrom^2'), symmetry=1, barrier=(22.149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.68747,0.063083,-4.71676e-05,1.5435e-08,-1.51135e-12,10271.3,23.1626], Tmin=(100,'K'), Tmax=(1227.37,'K')), NASAPolynomial(coeffs=[16.4739,0.0219359,-9.46965e-06,1.79658e-09,-1.26138e-13,5620.3,-59.3927], Tmin=(1227.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.3418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1C=C(C=O)C1(23610)',
    structure = SMILES('[CH2]C1C=C(C=O)C1'),
    E0 = (186.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63059,0.0394999,1.25142e-05,-4.13223e-08,1.74378e-11,22586.2,23.5199], Tmin=(100,'K'), Tmax=(1018.29,'K')), NASAPolynomial(coeffs=[12.598,0.0257541,-1.0451e-05,2.00441e-09,-1.45203e-13,18831.7,-37.0582], Tmin=(1018.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(Isobutyl)"""),
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
    label = '[C]=CC([CH2])[CH2](6428)',
    structure = SMILES('[C]=CC([CH2])[CH2]'),
    E0 = (919.993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.106372,'amu*angstrom^2'), symmetry=1, barrier=(2.4457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000998043,'amu*angstrom^2'), symmetry=1, barrier=(2.43852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.91323,'amu*angstrom^2'), symmetry=1, barrier=(66.9808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04587,0.0415728,-3.16589e-05,9.78588e-09,7.12838e-13,110721,22.3447], Tmin=(100,'K'), Tmax=(835.384,'K')), NASAPolynomial(coeffs=[8.28465,0.0196222,-6.46946e-06,1.03555e-09,-6.56747e-14,109403,-8.28617], Tmin=(835.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(919.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(CdCdJ2_triplet)"""),
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
    E0 = (435.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (478.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (504.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (693.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (564.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (482.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (576.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (685.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (607.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (826.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (510.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (526.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (754.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (598.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (884.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (813.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (622.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (565.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (560.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (551.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (458.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (460.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (460.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (720.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (740.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (673.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (637.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (750.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (592.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (443.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (723.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1056.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (890.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (618.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (585.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (672.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (597.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (611.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (523.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (524.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (443.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (987.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C=C(87)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C1CC1[C]=C[O](23586)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.264e+11,'s^-1'), n=0.35, Ea=(43.0952,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C(=C)C=C=C[O](22603)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(170.395,'m^3/(mol*s)'), n=1.5621, Ea=(11.2886,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2]C([CH2])C=C=C=O(23587)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(28)', 'C=CC=C=C[O](22346)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.0532,'m^3/(mol*s)'), n=2.065, Ea=(15.9938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;YJ] for rate rule [Cds-CdH_Cds-HH;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C=C(87)', '[CH]=C=C[O](8556)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.04,'cm^3/(mol*s)'), n=3.05, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CdsJ=Cdd]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C(C)=C[C]=C[O](22744)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([CH2])[CH]C#CO(23588)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C(C)[C]=C=C[O](22745)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.48e+08,'s^-1'), n=1.14, Ea=(172.799,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 205 used for R3H_SS_Cs;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH2])C#C[CH]O(23589)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C(C)[CH][C]=C=O(22746)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.42353e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C([CH2])=C[C]=CO(23590)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;O_rad_out;XH_out] for rate rule [R5H_SMMS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C=C[O](8556)', '[CH2][CH][CH2](6136)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(8)', '[CH2]C([CH2])=C[C]=C[O](23591)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH2]C([CH2])[C]=C=C[O](23592)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH2]C([CH2])[CH][C]=C=O(23593)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C([CH2])C=C1[CH]O1(23594)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C1[CH]C(=C[O])C1(23595)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.43734e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C([CH2])C1[C]=CO1(23596)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C1C=[C]C([O])C1(23597)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.98919e+10,'s^-1'), n=0.314866, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['C=C(C)C=C=C[O](22735)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C(=C)C=C=CO(23598)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C(C)C=C=C=O(22736)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][C]([CH2])C[C]=C[O](23599)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH2])[C]C=C[O](23600)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])[CH]C=[C][O](23601)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH2])=C[CH][CH][O](23602)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH2])=C[C]C[O](23603)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2][CH]CC=C=C[O](22491)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[O]C=C=CC1CC1(23604)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(T)(28)', '[CH2]C=C[C]=C[O](22677)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(T)(63)', 'C#C[CH]C([CH2])[CH2](19176)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(8)', '[CH]C([CH2])C=C=C[O](23605)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(8)', '[CH2]C([CH2])C#CC=O(23606)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH][CH2](6136)', 'C#CC=O(21959)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([CH2])[C]=CC=O(23607)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C([CH2])[CH]C=C=O(23608)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C([CH2])=CC=C[O](18261)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.2826e+06,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C1[CH][C]=COC1(23609)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.17824e+09,'s^-1'), n=0.529986, Ea=(88.0823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C(=C)C=CC=O(18274)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([CH2])C=C=C[O](22494)'],
    products = ['[CH2]C1C=C(C=O)C1(23610)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=O(373)', '[C]=CC([CH2])[CH2](6428)'],
    products = ['[CH2]C([CH2])C=C=C[O](22494)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4795',
    isomers = [
        '[CH2]C([CH2])C=C=C[O](22494)',
    ],
    reactants = [
        ('[CH2]C=C(87)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4795',
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

