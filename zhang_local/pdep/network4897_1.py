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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.362348,0.0800131,-8.91188e-05,5.1578e-08,-1.18964e-11,18624.5,29.091], Tmin=(100,'K'), Tmax=(1053.85,'K')), NASAPolynomial(coeffs=[14.8601,0.0249841,-1.07915e-05,2.02689e-09,-1.41383e-13,15568.8,-41.6122], Tmin=(1053.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CJC(C)C=O)"""),
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
    label = '[O]C=C=CC1CC1[O](25369)',
    structure = SMILES('[O]C=C=CC1CC1[O]'),
    E0 = (245.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796628,0.0492375,2.15034e-05,-7.68044e-08,3.66506e-11,29657,26.7306], Tmin=(100,'K'), Tmax=(939.345,'K')), NASAPolynomial(coeffs=[22.607,0.00997945,-1.4254e-06,2.33058e-10,-2.44476e-14,23194,-89.7183], Tmin=(939.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C1C=C=COC1[O](25370)',
    structure = SMILES('[CH2]C1C=C=COC1[O]'),
    E0 = (222.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20494,0.116577,-0.00015326,1.01834e-07,-2.65076e-11,26938.7,15.1875], Tmin=(100,'K'), Tmax=(944.742,'K')), NASAPolynomial(coeffs=[19.788,0.0276924,-1.21343e-05,2.24596e-09,-1.54161e-13,22972.2,-84.8977], Tmin=(944.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(Isobutyl) + radical(CCOJ)"""),
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
    label = 'C=C(C=O)C=C=C[O](25371)',
    structure = SMILES('C=C(C=O)C=C=C[O]'),
    E0 = (53.2164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25217,'amu*angstrom^2'), symmetry=1, barrier=(28.7898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25332,'amu*angstrom^2'), symmetry=1, barrier=(28.8164,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.147094,0.0704151,-4.65819e-05,-3.88273e-09,9.87133e-12,6552.05,25.1913], Tmin=(100,'K'), Tmax=(975.425,'K')), NASAPolynomial(coeffs=[22.9031,0.00963592,-3.15311e-06,6.33744e-10,-5.1253e-14,564.761,-91.9618], Tmin=(975.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.2164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C=O)C=C=C=O(25372)',
    structure = SMILES('[CH2]C(C=O)C=C=C=O'),
    E0 = (197.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,540,610,2055,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,4000],'cm^-1')),
        HinderedRotor(inertia=(0.230476,'amu*angstrom^2'), symmetry=1, barrier=(5.2991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000464094,'amu*angstrom^2'), symmetry=1, barrier=(5.26935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000460209,'amu*angstrom^2'), symmetry=1, barrier=(5.22524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62954,0.0565348,-7.0892e-05,5.36605e-08,-1.71091e-11,23878.3,11.749], Tmin=(100,'K'), Tmax=(755.124,'K')), NASAPolynomial(coeffs=[7.13203,0.0273842,-1.29802e-05,2.52724e-09,-1.78517e-13,23047.4,-13.2512], Tmin=(755.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds) + radical(CJC(C)C=O)"""),
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
    label = 'CC(C=O)=C[C]=C[O](25373)',
    structure = SMILES('CC(C=O)=C[C]=C[O]'),
    E0 = (58.7552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.231513,0.0765174,-7.48913e-05,3.67275e-08,-7.09628e-12,7207.64,26.821], Tmin=(100,'K'), Tmax=(1257.25,'K')), NASAPolynomial(coeffs=[17.6852,0.0209877,-8.63996e-06,1.59717e-09,-1.1073e-13,2818.92,-61.3784], Tmin=(1257.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.7552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(C=O)C=C=[C]O(25374)',
    structure = SMILES('[CH2]C(C=O)C=C=[C]O'),
    E0 = (252.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,3000,3100,440,815,1455,1000,1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.678683,'amu*angstrom^2'), symmetry=1, barrier=(15.6042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.678796,'amu*angstrom^2'), symmetry=1, barrier=(15.6069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.678794,'amu*angstrom^2'), symmetry=1, barrier=(15.6068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.678966,'amu*angstrom^2'), symmetry=1, barrier=(15.6108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0870892,0.0895812,-0.000122242,8.77623e-08,-2.49269e-11,30451.6,31.7995], Tmin=(100,'K'), Tmax=(864.56,'K')), NASAPolynomial(coeffs=[13.9365,0.0255055,-1.10717e-05,2.03892e-09,-1.38868e-13,28056.9,-33.0002], Tmin=(864.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CJC(C)C=O) + radical(C=CJO)"""),
)

species(
    label = 'CC([C]=C=C[O])C=O(25375)',
    structure = SMILES('CC([C]=C=C[O])C=O'),
    E0 = (181.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61419,0.0746273,-7.83552e-05,4.35335e-08,-9.75555e-12,21902.4,28.309], Tmin=(100,'K'), Tmax=(1076.64,'K')), NASAPolynomial(coeffs=[13.518,0.0266859,-1.15617e-05,2.17403e-09,-1.51664e-13,19123.9,-34.8973], Tmin=(1076.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'CC([C]=O)C=C=C[O](25376)',
    structure = SMILES('CC([C]=O)C=C=C[O]'),
    E0 = (101.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.349266,0.0728972,-7.01515e-05,3.41829e-08,-6.55836e-12,12398.3,29.9408], Tmin=(100,'K'), Tmax=(1268.58,'K')), NASAPolynomial(coeffs=[17.091,0.0201082,-7.73211e-06,1.37992e-09,-9.38239e-14,8150.66,-54.8109], Tmin=(1268.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C(C#C[CH]O)C=O(25377)',
    structure = SMILES('[CH2]C(C#C[CH]O)C=O'),
    E0 = (231.445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2100,2250,500,550,1380,1390,370,380,2900,435,215.929,2050.68],'cm^-1')),
        HinderedRotor(inertia=(2.00738,'amu*angstrom^2'), symmetry=1, barrier=(66.4193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243662,'amu*angstrom^2'), symmetry=1, barrier=(8.06212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243626,'amu*angstrom^2'), symmetry=1, barrier=(8.0621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243621,'amu*angstrom^2'), symmetry=1, barrier=(8.06199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00749,'amu*angstrom^2'), symmetry=1, barrier=(66.4191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.160536,0.0924753,-0.00014646,1.23118e-07,-3.94798e-11,27967.1,32.0225], Tmin=(100,'K'), Tmax=(919.976,'K')), NASAPolynomial(coeffs=[9.84641,0.0302523,-1.22186e-05,2.0793e-09,-1.31271e-13,27036,-9.27331], Tmin=(919.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CtOsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCsJOH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CC(C=O)C=[C][C]=O(25378)',
    structure = SMILES('CC(C=O)C=[C][C]=O'),
    E0 = (151.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27981,0.0741872,-4.96443e-05,-8.57334e-08,1.22985e-10,18255.7,26.8783], Tmin=(100,'K'), Tmax=(452.99,'K')), NASAPolynomial(coeffs=[8.51929,0.0376999,-1.96814e-05,3.88743e-09,-2.7302e-13,17318.3,-5.42254], Tmin=(452.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C(C=C=CO)=C[O](25379)',
    structure = SMILES('[CH2]C(C=C=CO)=C[O]'),
    E0 = (72.2858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.42989,0.0741445,-1.93668e-05,-5.92742e-08,3.73878e-11,8875.35,26.7673], Tmin=(100,'K'), Tmax=(904.884,'K')), NASAPolynomial(coeffs=[31.7552,-0.00493166,6.95647e-06,-1.48733e-09,9.9297e-14,463.293,-139.586], Tmin=(904.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.2858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)C=C=CO(25380)',
    structure = SMILES('[CH2]C([C]=O)C=C=CO'),
    E0 = (170.989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,540,610,2055,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.861565,'amu*angstrom^2'), symmetry=1, barrier=(19.8091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861825,'amu*angstrom^2'), symmetry=1, barrier=(19.8151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86221,'amu*angstrom^2'), symmetry=1, barrier=(19.8239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862068,'amu*angstrom^2'), symmetry=1, barrier=(19.8206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.283337,0.0891797,-0.000107168,6.34485e-08,-1.44685e-11,20723.8,31.1486], Tmin=(100,'K'), Tmax=(1085.34,'K')), NASAPolynomial(coeffs=[19.7133,0.0154838,-5.31788e-06,8.88157e-10,-5.83534e-14,16383.1,-66.9613], Tmin=(1085.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)CJ=O) + radical(CJC(C)C=O)"""),
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
    label = '[CH2]C(C=C=C[O])=C[O](25381)',
    structure = SMILES('[CH2]C(C=C=C[O])=C[O]'),
    E0 = (213.748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,540,610,2055,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49175,'amu*angstrom^2'), symmetry=1, barrier=(34.2984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49152,'amu*angstrom^2'), symmetry=1, barrier=(34.2929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0132216,0.0671327,-1.5933e-05,-5.28754e-08,3.24399e-11,25872.2,26.6808], Tmin=(100,'K'), Tmax=(915.722,'K')), NASAPolynomial(coeffs=[28.3918,-0.00105955,4.2272e-06,-9.07718e-10,5.76484e-14,18326.9,-120.65], Tmin=(915.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([C]=C=C[O])C=O(25382)',
    structure = SMILES('[CH2]C([C]=C=C[O])C=O'),
    E0 = (391.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.808729,'amu*angstrom^2'), symmetry=1, barrier=(18.5943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.804957,'amu*angstrom^2'), symmetry=1, barrier=(18.5075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.808165,'amu*angstrom^2'), symmetry=1, barrier=(18.5813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331291,0.0845215,-0.000113308,7.97758e-08,-2.23224e-11,47227.9,29.6811], Tmin=(100,'K'), Tmax=(875.104,'K')), NASAPolynomial(coeffs=[13.3347,0.0250846,-1.14283e-05,2.16307e-09,-1.50091e-13,44952,-31.3181], Tmin=(875.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([C]=O)C=C=C[O](25383)',
    structure = SMILES('[CH2]C([C]=O)C=C=C[O]'),
    E0 = (312.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.856626,'amu*angstrom^2'), symmetry=1, barrier=(19.6955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.857978,'amu*angstrom^2'), symmetry=1, barrier=(19.7266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.856658,'amu*angstrom^2'), symmetry=1, barrier=(19.6963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.364241,0.0793522,-9.34409e-05,5.5821e-08,-1.31103e-11,37710.8,30.2399], Tmin=(100,'K'), Tmax=(1042.57,'K')), NASAPolynomial(coeffs=[15.7205,0.0204348,-8.67279e-06,1.61612e-09,-1.12341e-13,34508.8,-44.4855], Tmin=(1042.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)CJ=O) + radical(C=COJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=O)C=[C][C]=O(25384)',
    structure = SMILES('[CH2]C(C=O)C=[C][C]=O'),
    E0 = (361.609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3000,3100,440,815,1455,1000,1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.608667,'amu*angstrom^2'), symmetry=1, barrier=(13.9945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60636,'amu*angstrom^2'), symmetry=1, barrier=(13.9414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.605842,'amu*angstrom^2'), symmetry=1, barrier=(13.9295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.606872,'amu*angstrom^2'), symmetry=1, barrier=(13.9532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0515362,0.10156,-0.000182,1.69565e-07,-6.07835e-11,43625.3,31.7449], Tmin=(100,'K'), Tmax=(828.32,'K')), NASAPolynomial(coeffs=[8.68967,0.0355074,-1.9212e-05,3.79801e-09,-2.65001e-13,42995.1,-3.84288], Tmin=(828.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=O) + radical(CJC(C)C=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C(C=O)C=C1[CH]O1(25385)',
    structure = SMILES('[CH2]C(C=O)C=C1[CH]O1'),
    E0 = (179.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.490302,0.0623574,-2.46385e-05,-2.37803e-08,1.64869e-11,21782.2,27.207], Tmin=(100,'K'), Tmax=(957.769,'K')), NASAPolynomial(coeffs=[20.7164,0.0135378,-4.01666e-06,7.31209e-10,-5.60534e-14,16272.6,-78.0362], Tmin=(957.769,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(methyleneoxirane) + radical(CJC(C)C=O) + radical(C=CCJO)"""),
)

species(
    label = '[O]C=C1[CH]C(C=O)C1(25386)',
    structure = SMILES('[O]C=C1[CH]C(C=O)C1'),
    E0 = (120.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.339,-0.0106509,0.000122985,-1.21401e-07,2.82317e-11,14162.1,-12.4098], Tmin=(100,'K'), Tmax=(1711.79,'K')), NASAPolynomial(coeffs=[76.1395,0.0219783,-6.89343e-05,1.69521e-08,-1.26403e-12,-35673,-444.994], Tmin=(1711.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsOsH) + ring(methylenecyclobutane) + radical(C=COJ) + radical(CCJCC=O)"""),
)

species(
    label = '[O]C=C=CC1[CH]OC1(25387)',
    structure = SMILES('[O]C=C=CC1[CH]OC1'),
    E0 = (229.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.718639,0.0540884,1.12912e-05,-7.50759e-08,4.07318e-11,27754.2,25.5498], Tmin=(100,'K'), Tmax=(874.239,'K')), NASAPolynomial(coeffs=[22.5068,0.00648341,3.60514e-06,-1.06744e-09,8.05237e-14,21954.2,-88.0206], Tmin=(874.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Oxetane) + radical(C=COJ) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C(C=O)C1[C]=CO1(25388)',
    structure = SMILES('[CH2]C(C=O)C1[C]=CO1'),
    E0 = (274.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.327852,0.0640731,-2.23273e-05,-3.09224e-08,2.01143e-11,33190.4,28.1649], Tmin=(100,'K'), Tmax=(950.593,'K')), NASAPolynomial(coeffs=[22.7336,0.0106734,-2.57403e-06,4.65896e-10,-3.88262e-14,27083.6,-88.5102], Tmin=(950.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(CJC(C)C=O)"""),
)

species(
    label = '[O]C1[C]=CC(C=O)C1(25389)',
    structure = SMILES('[O]C1[C]=CC(C=O)C1'),
    E0 = (223.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32417,0.0476189,-7.64945e-06,-2.15865e-08,1.0607e-11,26996.8,25.9062], Tmin=(100,'K'), Tmax=(1053.27,'K')), NASAPolynomial(coeffs=[13.3816,0.025845,-1.08434e-05,2.08397e-09,-1.49807e-13,23124.7,-39.2136], Tmin=(1053.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH2]C1[CH]OOC=C=C1(25390)',
    structure = SMILES('[CH2]C1[CH]OOC=C=C1'),
    E0 = (540.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0693,0.0539737,-1.81841e-05,-1.33946e-08,8.42998e-12,65121.5,22.4124], Tmin=(100,'K'), Tmax=(1036.18,'K')), NASAPolynomial(coeffs=[13.5674,0.0273255,-1.08742e-05,2.01868e-09,-1.42251e-13,61371.9,-43.9231], Tmin=(1036.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(540.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCsJOOC) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(C=O)C=C=CO(25391)',
    structure = SMILES('C=C(C=O)C=C=CO'),
    E0 = (-88.2462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.261088,0.0773556,-4.99312e-05,-1.0087e-08,1.45694e-11,-10445.2,25.2458], Tmin=(100,'K'), Tmax=(945.532,'K')), NASAPolynomial(coeffs=[26.1011,0.00603994,-5.80957e-07,9.08862e-11,-1.26297e-14,-17227.8,-109.964], Tmin=(945.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.2462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC(C=O)C=C=C=O(25392)',
    structure = SMILES('CC(C=O)C=C=C=O'),
    E0 = (-12.6515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04231,0.0452261,-3.13019e-05,1.13868e-08,-1.77143e-12,-1453.39,9.90043], Tmin=(100,'K'), Tmax=(1411.73,'K')), NASAPolynomial(coeffs=[8.44839,0.0270751,-1.20162e-05,2.27946e-09,-1.58646e-13,-3262.13,-23.2142], Tmin=(1411.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.6515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C(=C[O])C[C]=C[O](25393)',
    structure = SMILES('[CH2]C(=C[O])C[C]=C[O]'),
    E0 = (304.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17297,'amu*angstrom^2'), symmetry=1, barrier=(26.9689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1736,'amu*angstrom^2'), symmetry=1, barrier=(26.9834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17328,'amu*angstrom^2'), symmetry=1, barrier=(26.9761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.132994,0.0674605,-2.41013e-05,-3.49708e-08,2.33238e-11,36763.4,30.7534], Tmin=(100,'K'), Tmax=(930.831,'K')), NASAPolynomial(coeffs=[24.4587,0.00787601,-5.16828e-07,1.53078e-11,-5.72556e-15,30287.5,-95.3206], Tmin=(930.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=COJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=C[C]=C[O])C[O](25394)',
    structure = SMILES('[CH2]C(=C[C]=C[O])C[O]'),
    E0 = (341.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,402.631,403.532,403.574,404.345],'cm^-1')),
        HinderedRotor(inertia=(0.645674,'amu*angstrom^2'), symmetry=1, barrier=(75.1521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.651054,'amu*angstrom^2'), symmetry=1, barrier=(75.1595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.648825,'amu*angstrom^2'), symmetry=1, barrier=(75.296,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424019,0.0725963,-7.00777e-05,3.60259e-08,-7.38068e-12,41162.5,30.4323], Tmin=(100,'K'), Tmax=(1187.11,'K')), NASAPolynomial(coeffs=[14.8981,0.0238257,-8.45285e-06,1.41834e-09,-9.25088e-14,37726,-41.8794], Tmin=(1187.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([C]C=C[O])C=O(25395)',
    structure = SMILES('[CH2]C([C]C=C[O])C=O'),
    E0 = (438.754,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,332.669,332.669,332.67,332.67],'cm^-1')),
        HinderedRotor(inertia=(0.222186,'amu*angstrom^2'), symmetry=1, barrier=(17.4491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308016,'amu*angstrom^2'), symmetry=1, barrier=(24.1896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222186,'amu*angstrom^2'), symmetry=1, barrier=(17.4491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222184,'amu*angstrom^2'), symmetry=1, barrier=(17.4491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.098952,0.0817703,-8.88294e-05,4.86299e-08,-1.04423e-11,52913.8,31.9905], Tmin=(100,'K'), Tmax=(1139.28,'K')), NASAPolynomial(coeffs=[17.5143,0.0206255,-8.325e-06,1.52164e-09,-1.04977e-13,48945.6,-54.2992], Tmin=(1139.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CCJ2_triplet) + radical(CJC(C)C=O) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH][CH][C]=O)C=O(25396)',
    structure = SMILES('[CH2]C([CH][CH][C]=O)C=O'),
    E0 = (374.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3000,3050,390,425,1340,1360,335,370,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,360.65,2512.78],'cm^-1')),
        HinderedRotor(inertia=(0.0819972,'amu*angstrom^2'), symmetry=1, barrier=(7.56907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00129665,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0012969,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.081995,'amu*angstrom^2'), symmetry=1, barrier=(7.5686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712935,'amu*angstrom^2'), symmetry=1, barrier=(65.8422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576706,0.0797189,-0.000107247,8.21718e-08,-2.55883e-11,45182.6,33.7376], Tmin=(100,'K'), Tmax=(794.774,'K')), NASAPolynomial(coeffs=[10.1303,0.0308987,-1.37138e-05,2.54576e-09,-1.73825e-13,43687.3,-10.0116], Tmin=(794.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CJC(C)C=O) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=C=C[O])C[O](25397)',
    structure = SMILES('[CH2]C([C]=C=C[O])C[O]'),
    E0 = (533.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,335.528,335.555,335.564,335.565],'cm^-1')),
        HinderedRotor(inertia=(0.00149769,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123495,'amu*angstrom^2'), symmetry=1, barrier=(9.86805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.930216,'amu*angstrom^2'), symmetry=1, barrier=(74.3251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.531022,0.0779739,-9.11309e-05,5.84432e-08,-1.51124e-11,64322,32.4842], Tmin=(100,'K'), Tmax=(940.79,'K')), NASAPolynomial(coeffs=[12.2569,0.0281182,-1.16403e-05,2.11398e-09,-1.43755e-13,62115.7,-23.3707], Tmin=(940.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(Cds_S) + radical(Isobutyl) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([C]=O)C[C]=C[O](25398)',
    structure = SMILES('[CH2]C([C]=O)C[C]=C[O]'),
    E0 = (387.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,394.187,394.187,394.187],'cm^-1')),
        HinderedRotor(inertia=(0.125018,'amu*angstrom^2'), symmetry=1, barrier=(13.7849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125018,'amu*angstrom^2'), symmetry=1, barrier=(13.7849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125018,'amu*angstrom^2'), symmetry=1, barrier=(13.7849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125018,'amu*angstrom^2'), symmetry=1, barrier=(13.7849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282675,0.0806978,-9.14217e-05,5.32073e-08,-1.22366e-11,46688.5,33.0583], Tmin=(100,'K'), Tmax=(1061.8,'K')), NASAPolynomial(coeffs=[15.683,0.0226825,-9.46466e-06,1.74997e-09,-1.21122e-13,43418,-42.163], Tmin=(1061.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CJC(C)C=O) + radical(C=COJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C([CH]C=C[O])=C[O](25399)',
    structure = SMILES('[CH2]C([CH]C=C[O])=C[O]'),
    E0 = (168.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,299.388,299.728,299.825,300.297],'cm^-1')),
        HinderedRotor(inertia=(0.593617,'amu*angstrom^2'), symmetry=1, barrier=(37.9514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.597153,'amu*angstrom^2'), symmetry=1, barrier=(37.9528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.595758,'amu*angstrom^2'), symmetry=1, barrier=(37.9444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483379,0.0504761,3.83671e-05,-1.09151e-07,5.22697e-11,20388.5,28.6933], Tmin=(100,'K'), Tmax=(917.438,'K')), NASAPolynomial(coeffs=[27.8682,0.00124827,4.12779e-06,-9.03818e-10,5.51184e-14,12410.7,-117.157], Tmin=(917.438,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CCJC=C) + radical(C=COJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=C=C[O])[CH]O(25400)',
    structure = SMILES('[CH2]C([C]=C=C[O])[CH]O'),
    E0 = (488.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,540,610,2055,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,330.076,330.338,330.545],'cm^-1')),
        HinderedRotor(inertia=(1.4431,'amu*angstrom^2'), symmetry=1, barrier=(111.654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179931,'amu*angstrom^2'), symmetry=1, barrier=(13.9597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180771,'amu*angstrom^2'), symmetry=1, barrier=(13.9639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180301,'amu*angstrom^2'), symmetry=1, barrier=(13.9615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.183705,0.0905117,-0.000114055,7.13459e-08,-1.68537e-11,58889.3,33.7332], Tmin=(100,'K'), Tmax=(885.401,'K')), NASAPolynomial(coeffs=[17.6966,0.0182728,-6.13917e-06,9.83019e-10,-6.19329e-14,55388.4,-52.243], Tmin=(885.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([C]=O)[CH]C=C[O](25401)',
    structure = SMILES('[CH2]C([C]=O)[CH]C=C[O]'),
    E0 = (349.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,428.286,432.992,433.128],'cm^-1')),
        HinderedRotor(inertia=(0.000902929,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11859,'amu*angstrom^2'), symmetry=1, barrier=(15.8667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117489,'amu*angstrom^2'), symmetry=1, barrier=(15.8492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119459,'amu*angstrom^2'), symmetry=1, barrier=(15.8706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0547893,0.0762222,-7.64858e-05,3.80792e-08,-7.34697e-12,42141.5,35.8999], Tmin=(100,'K'), Tmax=(1275.5,'K')), NASAPolynomial(coeffs=[19.3447,0.0157292,-5.34601e-06,8.96819e-10,-5.92269e-14,37220.6,-61.8568], Tmin=(1275.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CC(C)CJ=O) + radical(C=COJ) + radical(CJC(C)C=O) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C(C=[C]C[O])=C[O](25402)',
    structure = SMILES('[CH2]C(C=[C]C[O])=C[O]'),
    E0 = (414.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,395.315,395.406,395.504,395.546],'cm^-1')),
        HinderedRotor(inertia=(0.14677,'amu*angstrom^2'), symmetry=1, barrier=(16.2793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146663,'amu*angstrom^2'), symmetry=1, barrier=(16.2771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276068,'amu*angstrom^2'), symmetry=1, barrier=(30.6475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296051,0.075771,-7.5836e-05,3.87854e-08,-7.82984e-12,50037.4,30.3368], Tmin=(100,'K'), Tmax=(1206.28,'K')), NASAPolynomial(coeffs=[16.6579,0.0215148,-8.36813e-06,1.49805e-09,-1.02003e-13,46090,-51.6682], Tmin=(1206.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(CCOJ) + radical(Cds_S) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)C=[C]C[O](25403)',
    structure = SMILES('[CH2]C([C]=O)C=[C]C[O]'),
    E0 = (513.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180,180,2328.21],'cm^-1')),
        HinderedRotor(inertia=(0.217115,'amu*angstrom^2'), symmetry=1, barrier=(4.99189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217134,'amu*angstrom^2'), symmetry=1, barrier=(4.99235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217128,'amu*angstrom^2'), symmetry=1, barrier=(4.9922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217094,'amu*angstrom^2'), symmetry=1, barrier=(4.99142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.603278,0.0889456,-0.000157466,1.54373e-07,-5.72041e-11,61879,34.1411], Tmin=(100,'K'), Tmax=(856.115,'K')), NASAPolynomial(coeffs=[3.09629,0.0444964,-2.21158e-05,4.22051e-09,-2.88372e-13,62654.2,29.5215], Tmin=(856.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCOJ) + radical(Cds_S) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C([CH][C]=C=O)C[O](25404)',
    structure = SMILES('[CH2]C([CH][C]=C=O)C[O]'),
    E0 = (462.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,266.915,267.239,267.307,2174.77],'cm^-1')),
        HinderedRotor(inertia=(1.09076,'amu*angstrom^2'), symmetry=1, barrier=(55.2015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182707,'amu*angstrom^2'), symmetry=1, barrier=(9.26035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164468,'amu*angstrom^2'), symmetry=1, barrier=(55.1991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164487,'amu*angstrom^2'), symmetry=1, barrier=(55.2094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619138,0.0829256,-0.000125255,1.11342e-07,-3.91332e-11,55710.2,32.9495], Tmin=(100,'K'), Tmax=(845.747,'K')), NASAPolynomial(coeffs=[6.17846,0.0396744,-1.84695e-05,3.45885e-09,-2.35378e-13,55376.4,10.6459], Tmin=(845.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(Isobutyl) + radical(CCJC(C)=C=O) + radical(CCCJ=C=O) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([CH]O)[CH][C]=C=O(25405)',
    structure = SMILES('[CH2]C([CH]O)[CH][C]=C=O'),
    E0 = (416.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0134918,0.0940862,-0.000142892,1.16597e-07,-3.72008e-11,50273,33.8129], Tmin=(100,'K'), Tmax=(873.477,'K')), NASAPolynomial(coeffs=[11.4469,0.0301307,-1.31466e-05,2.37073e-09,-1.57153e-13,48718,-17.2674], Tmin=(873.477,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCJC(C)=C=O) + radical(CCsJOH) + radical(Isobutyl) + radical(CCCJ=C=O)"""),
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
    label = '[CH2]CC=C=C[O](22300)',
    structure = SMILES('[CH2]CC=C=C[O]'),
    E0 = (261.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.978739,'amu*angstrom^2'), symmetry=1, barrier=(22.5031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.980603,'amu*angstrom^2'), symmetry=1, barrier=(22.546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3802.23,'J/mol'), sigma=(6.27918,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.90 K, Pc=34.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44142,0.0451849,-1.09233e-05,-2.37761e-08,1.39002e-11,31610.3,23.1768], Tmin=(100,'K'), Tmax=(959.35,'K')), NASAPolynomial(coeffs=[15.1908,0.0145156,-4.65246e-06,8.31917e-10,-6.07325e-14,27745.4,-48.9793], Tmin=(959.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(RCCJ)"""),
)

species(
    label = '[O]C=C=CCC=C[O](22588)',
    structure = SMILES('[O]C=C=CCC=C[O]'),
    E0 = (94.6684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371633,0.0609434,-8.25709e-06,-4.93002e-08,2.78157e-11,11534,29.2454], Tmin=(100,'K'), Tmax=(935.092,'K')), NASAPolynomial(coeffs=[23.8145,0.00848651,-8.25121e-07,9.47677e-11,-1.27022e-14,5058.87,-93.4597], Tmin=(935.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.6684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = 'C=C[CH]OC=C=C[O](22592)',
    structure = SMILES('C=C[CH]OC=C=C[O]'),
    E0 = (108.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39351,'amu*angstrom^2'), symmetry=1, barrier=(32.0395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4005,'amu*angstrom^2'), symmetry=1, barrier=(32.2003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3957,'amu*angstrom^2'), symmetry=1, barrier=(32.0899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194805,0.0660784,-2.13151e-05,-3.5764e-08,2.28569e-11,13171.4,28.5084], Tmin=(100,'K'), Tmax=(939.703,'K')), NASAPolynomial(coeffs=[23.8309,0.0096062,-1.62765e-06,2.5359e-10,-2.32835e-14,6780.43,-94.4212], Tmin=(939.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJ(O)C) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C=C=C[O])=CO(25407)',
    structure = SMILES('[CH2]C(C=C=C[O])=CO'),
    E0 = (72.2858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35624,'amu*angstrom^2'), symmetry=1, barrier=(31.1827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35316,'amu*angstrom^2'), symmetry=1, barrier=(31.1118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35546,'amu*angstrom^2'), symmetry=1, barrier=(31.1646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.42989,0.0741445,-1.93668e-05,-5.92742e-08,3.73878e-11,8875.35,26.7673], Tmin=(100,'K'), Tmax=(904.884,'K')), NASAPolynomial(coeffs=[31.7552,-0.00493166,6.95647e-06,-1.48733e-09,9.9297e-14,463.293,-139.586], Tmin=(904.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.2858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C1[CH]OC([CH]1)=C[O](25408)',
    structure = SMILES('[CH2]C1[CH]OC([CH]1)=C[O]'),
    E0 = (349.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.809148,0.0384042,7.4744e-05,-1.51875e-07,6.93042e-11,42120.4,29.5574], Tmin=(100,'K'), Tmax=(908.761,'K')), NASAPolynomial(coeffs=[29.7064,-0.00506356,8.29336e-06,-1.74451e-09,1.13009e-13,33411.1,-126.112], Tmin=(908.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(CCJCO) + radical(Isobutyl) + radical(CCsJOC(O)) + radical(C=COJ)"""),
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
    label = '[O]C=C=CC=C[O](23718)',
    structure = SMILES('[O]C=C=CC=C[O]'),
    E0 = (99.8399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.45304,'amu*angstrom^2'), symmetry=1, barrier=(33.4082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.811775,0.0498515,5.24762e-06,-6.54747e-08,3.53354e-11,12141.7,22.8384], Tmin=(100,'K'), Tmax=(912.976,'K')), NASAPolynomial(coeffs=[25.3204,-0.00384644,5.27499e-06,-1.09203e-09,7.01856e-14,5429.3,-105.423], Tmin=(912.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.8399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[CH]=C=CC([CH2])C=O(19966)',
    structure = SMILES('[CH]=C=CC([CH2])C=O'),
    E0 = (375.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.532983,'amu*angstrom^2'), symmetry=1, barrier=(12.2543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.531765,'amu*angstrom^2'), symmetry=1, barrier=(12.2263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.058,'amu*angstrom^2'), symmetry=1, barrier=(24.3256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835179,0.0742316,-0.000103291,8.24322e-08,-2.64427e-11,45280.1,26.7472], Tmin=(100,'K'), Tmax=(843.498,'K')), NASAPolynomial(coeffs=[8.98071,0.0294866,-1.2842e-05,2.34628e-09,-1.58011e-13,44123.6,-9.87394], Tmin=(843.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJC(C)C=O) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C(C=O)C=C=C[O](25409)',
    structure = SMILES('[CH]C(C=O)C=C=C[O]'),
    E0 = (391.467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,487.353,708.088,708.111],'cm^-1')),
        HinderedRotor(inertia=(0.605822,'amu*angstrom^2'), symmetry=1, barrier=(13.929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.91099,'amu*angstrom^2'), symmetry=1, barrier=(20.9455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.039141,'amu*angstrom^2'), symmetry=1, barrier=(13.9291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.43208,0.0750423,-8.00948e-05,4.32092e-08,-9.17802e-12,47214,29.4896], Tmin=(100,'K'), Tmax=(1148.51,'K')), NASAPolynomial(coeffs=[16.2165,0.0200688,-8.29735e-06,1.53355e-09,-1.06365e-13,43588.3,-48.8468], Tmin=(1148.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C1C=C(C=O)C1[O](25410)',
    structure = SMILES('[CH2]C1C=C(C=O)C1[O]'),
    E0 = (244.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09596,0.050088,-4.52485e-06,-3.13835e-08,1.54818e-11,29528.7,27.8206], Tmin=(100,'K'), Tmax=(1014.33,'K')), NASAPolynomial(coeffs=[16.0864,0.0220444,-9.00167e-06,1.75792e-09,-1.29543e-13,24889.3,-52.5921], Tmin=(1014.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(Isobutyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(C#CC=O)C=O(25411)',
    structure = SMILES('[CH2]C(C#CC=O)C=O'),
    E0 = (116.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,2100,2250,500,550,3000,3100,440,815,1455,1000,202.212],'cm^-1')),
        HinderedRotor(inertia=(4.12011,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292545,'amu*angstrom^2'), symmetry=1, barrier=(8.85802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303522,'amu*angstrom^2'), symmetry=1, barrier=(8.84451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.99036,'amu*angstrom^2'), symmetry=1, barrier=(88.8749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.691139,0.0789883,-0.000116008,9.73357e-08,-3.26791e-11,14125.3,28.2212], Tmin=(100,'K'), Tmax=(830.572,'K')), NASAPolynomial(coeffs=[8.47758,0.032249,-1.49103e-05,2.79441e-09,-1.90732e-13,13150.6,-5.97983], Tmin=(830.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([C]=CC=O)C=O(25412)',
    structure = SMILES('[CH2]C([C]=CC=O)C=O'),
    E0 = (194.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,182.55],'cm^-1')),
        HinderedRotor(inertia=(0.395726,'amu*angstrom^2'), symmetry=1, barrier=(9.49828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402541,'amu*angstrom^2'), symmetry=1, barrier=(9.47063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387146,'amu*angstrom^2'), symmetry=1, barrier=(9.45558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.409307,'amu*angstrom^2'), symmetry=1, barrier=(9.52207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.446376,0.089238,-0.000147973,1.40339e-07,-5.23031e-11,23554.4,30.2007], Tmin=(100,'K'), Tmax=(809.238,'K')), NASAPolynomial(coeffs=[5.60341,0.0425648,-2.21959e-05,4.37496e-09,-3.06753e-13,23413.3,10.698], Tmin=(809.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=O)C=C[C]=O(25413)',
    structure = SMILES('[CH2]C(C=O)C=C[C]=O'),
    E0 = (117.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.109808,0.0963975,-0.000161232,1.48189e-07,-5.33565e-11,14279.3,30.4614], Tmin=(100,'K'), Tmax=(816.261,'K')), NASAPolynomial(coeffs=[7.84971,0.0390837,-2.02866e-05,3.97908e-09,-2.77779e-13,13661.6,-1.35195], Tmin=(816.261,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=O)=CC=C[O](25414)',
    structure = SMILES('[CH2]C(C=O)=CC=C[O]'),
    E0 = (16.1113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.226829,0.0686073,-3.90368e-05,-9.08298e-09,1.08372e-11,2086.42,26.7384], Tmin=(100,'K'), Tmax=(986.381,'K')), NASAPolynomial(coeffs=[21.5288,0.0142942,-5.21309e-06,1.01944e-09,-7.78053e-14,-3676.13,-83.6478], Tmin=(986.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.1113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C=O)CJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([C]=O)C=CC=O(25415)',
    structure = SMILES('[CH2]C([C]=O)C=CC=O'),
    E0 = (115.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597235,0.0825487,-0.000122254,1.0824e-07,-3.94974e-11,14032.4,30.3456], Tmin=(100,'K'), Tmax=(752.659,'K')), NASAPolynomial(coeffs=[7.28697,0.0391647,-2.0185e-05,4.00819e-09,-2.84221e-13,13247.2,1.44601], Tmin=(752.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)CJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C1[CH]OC(C=O)=C1(25416)',
    structure = SMILES('[CH2]C1[CH]OC(C=O)=C1'),
    E0 = (117.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71675,0.0575312,-1.24975e-05,-3.85374e-08,2.33178e-11,14261.6,22.5713], Tmin=(100,'K'), Tmax=(913.915,'K')), NASAPolynomial(coeffs=[20.0188,0.0119771,-1.62012e-06,1.33358e-10,-9.36869e-15,9107.88,-77.7064], Tmin=(913.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + ring(2,3-Dihydrofuran) + radical(CCsJOC(O)) + radical(Isobutyl)"""),
)

species(
    label = 'O=CC1C=[C][CH]OC1(25417)',
    structure = SMILES('O=CC1C=[C][CH]OC1'),
    E0 = (115.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50886,0.0504627,-2.85663e-05,8.32264e-09,-1.01383e-12,13931.5,22.5944], Tmin=(100,'K'), Tmax=(1795.42,'K')), NASAPolynomial(coeffs=[10.7416,0.0298928,-1.13807e-05,1.94126e-09,-1.2525e-13,10616.2,-27.3515], Tmin=(1795.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(36dihydro2hpyran) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = 'C=C(C=O)C=CC=O(25418)',
    structure = SMILES('C=C(C=O)C=CC=O'),
    E0 = (-143.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71203,0.0694775,-5.93232e-05,2.41092e-08,-3.90413e-12,-17141.2,24.1149], Tmin=(100,'K'), Tmax=(1454.04,'K')), NASAPolynomial(coeffs=[17.3936,0.0235873,-1.19825e-05,2.40383e-09,-1.72239e-13,-21992.4,-62.6089], Tmin=(1454.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'O=CC1=CC(C=O)C1(25419)',
    structure = SMILES('O=CC1=CC(C=O)C1'),
    E0 = (-99.8132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24895,0.0537625,-3.18406e-05,7.78255e-09,-5.84595e-13,-11900.1,24.6169], Tmin=(100,'K'), Tmax=(1630.24,'K')), NASAPolynomial(coeffs=[16.6641,0.0234639,-1.08857e-05,2.04445e-09,-1.3881e-13,-17926,-60.3525], Tmin=(1630.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.8132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene)"""),
)

species(
    label = '[C]=CC([CH2])C=O(13854)',
    structure = SMILES('[C]=CC([CH2])C=O'),
    E0 = (638.616,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.259038,'amu*angstrom^2'), symmetry=1, barrier=(5.95578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256346,'amu*angstrom^2'), symmetry=1, barrier=(5.89389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.452631,'amu*angstrom^2'), symmetry=1, barrier=(10.4069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40526,0.0642923,-0.000102887,9.36695e-08,-3.33965e-11,76894.2,23.5902], Tmin=(100,'K'), Tmax=(837.419,'K')), NASAPolynomial(coeffs=[5.84108,0.028897,-1.40378e-05,2.67765e-09,-1.83957e-13,76649.5,5.95134], Tmin=(837.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(CJC(C)C=O)"""),
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
    E0 = (153.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (196.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (245.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (308.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (282.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (412.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (201.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (222.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (252.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (429.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (326.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (271.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (546.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (229.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (244.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (253.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (360.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (340.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (425.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (603.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (524.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (573.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (340.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (279.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (340.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (279.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (270.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (540.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (178.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (178.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (332.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (363.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (461.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (397.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (622.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (465.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (257.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (527.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (374.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (454.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (538.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (487.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (441.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (867.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (313.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (276.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (422.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (216.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (499.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (515.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (782.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (603.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (244.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (338.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (191.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (391.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (316.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (330.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (198.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (209.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (241.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (242.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (162.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (706.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['C=CC=O(5269)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[O]C=[C]C1CC1C=O(25368)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.32e+10,'s^-1'), n=0.35, Ea=(43.0952,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[O]C=C=CC1CC1[O](25369)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.93521e+09,'s^-1'), n=0.743095, Ea=(91.6988,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 87.6 to 91.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C1C=C=COC1[O](25370)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.04811e+11,'s^-1'), n=0.222, Ea=(155.084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;multiplebond_intra;radadd_intra_O] + [R7_SMMS;multiplebond_intra;radadd_intra] for rate rule [R7_SMMS;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C=C(C=O)C=C=C[O](25371)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(452.51,'m^3/(mol*s)'), n=1.51729, Ea=(17.7291,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-TwoDe_Cds;HJ] for rate rule [Cds-CdCO_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', '[CH2]C(C=O)C=C=C=O(25372)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=CC=O(5269)', '[CH]=C=C[O](8556)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00503627,'m^3/(mol*s)'), n=2.41968, Ea=(13.2583,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeH_Cds;CJ] for rate rule [Cds-COH_Cds;CdsJ=Cdd]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=O(373)', 'C=CC=C=C[O](22346)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00684716,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;CJ] for rate rule [Cds-CdH_Cds-HH;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['CC(C=O)=C[C]=C[O](25373)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C=O)C=C=[C]O(25374)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['CC([C]=C=C[O])C=O(25375)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.24e+08,'s^-1'), n=1.14, Ea=(172.799,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 205 used for R3H_SS_Cs;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['CC([C]=O)C=C=C[O](25376)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_2H;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C#C[CH]O)C=O(25377)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['CC(C=O)C=[C][C]=O(25378)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.21176e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C(C=C=CO)=C[O](25379)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;O_rad_out;XH_out] for rate rule [R5H_SMMS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([C]=O)C=C=CO(25380)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8024.52,'s^-1'), n=1.96653, Ea=(82.3883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Y_rad_out;XH_out] for rate rule [R6H;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C=C[O](5266)', '[CH]=C=C[O](8556)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=O(373)', '[CH2]C=C[C]=C[O](22677)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.7839e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2]C(C=C=C[O])=C[O](25381)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[CH2]C([C]=C=C[O])C=O(25382)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2]C([C]=O)C=C=C[O](25383)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH2]C(C=O)C=[C][C]=O(25384)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C(C=O)C=C1[CH]O1(25385)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[O]C=C1[CH]C(C=O)C1(25386)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[O]C=C=CC1[CH]OC1(25387)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C(C=O)C1[C]=CO1(25388)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[O]C1[C]=CC(C=O)C1(25389)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.4946e+10,'s^-1'), n=0.314866, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C1[CH]OOC=C=C1(25390)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(386.735,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 383.0 to 386.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['C=C(C=O)C=C=CO(25391)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['CC(C=O)C=C=C=O(25392)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(=C[O])C[C]=C[O](25393)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(=C[C]=C[O])C[O](25394)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([C]C=C[O])C=O(25395)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH][CH][C]=O)C=O(25396)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([C]=C=C[O])C[O](25397)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([C]=O)C[C]=C[O](25398)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH]C=C[O])=C[O](25399)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([C]=C=C[O])[CH]O(25400)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([C]=O)[CH]C=C[O](25401)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(C=[C]C[O])=C[O](25402)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([C]=O)C=[C]C[O](25403)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([CH][C]=C=O)C[O](25404)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([CH]O)[CH][C]=C=O(25405)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[C-]#[O+](374)', '[CH2]CC=C=C[O](22300)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[O]C=[C]C=CCC=O(25406)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C[CH]OC=C=C[O](22592)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C1[CH]OC([CH]1)=C[O](25408)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction50',
    reactants = ['CH2(T)(28)', '[O]C=C=CC=C[O](23718)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['O(T)(63)', '[CH]=C=CC([CH2])C=O(19966)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction52',
    reactants = ['H(8)', '[CH]C(C=O)C=C=C[O](25409)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C1C=C(C=O)C1[O](25410)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(90.7833,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonylbond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 86.4 to 90.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction54',
    reactants = ['H(8)', '[CH2]C(C#CC=O)C=O(25411)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C=C[O](5266)', 'C#CC=O(21959)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C([C]=CC=O)C=O(25412)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C(C=O)C=C[C]=O(25413)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C(C=O)=CC=C[O](25414)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(8.2826e+06,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C([C]=O)C=CC=O(25415)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleDe;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['[CH2]C1[CH]OC(C=O)=C1(25416)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(8.72e+09,'s^-1'), n=0.186, Ea=(55.7727,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleDe] for rate rule [R5_DS_CO;carbonyl_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['O=CC1C=[C][CH]OC1(25417)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(88.0823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['C=C(C=O)C=CC=O(25418)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C(C=O)C=C=C[O](22596)'],
    products = ['O=CC1=CC(C=O)C1(25419)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]=O(373)', '[C]=CC([CH2])C=O(13854)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4897',
    isomers = [
        '[CH2]C(C=O)C=C=C[O](22596)',
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
    label = '4897',
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

