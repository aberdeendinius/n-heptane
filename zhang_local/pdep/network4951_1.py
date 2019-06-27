species(
    label = '[CH]=C=COOC=C=[CH](22648)',
    structure = SMILES('[CH]=C=COOC=C=[CH]'),
    E0 = (653.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,540,563.333,586.667,610,1970,2140,350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(1.46479,'amu*angstrom^2'), symmetry=1, barrier=(33.6785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46525,'amu*angstrom^2'), symmetry=1, barrier=(33.6891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.25457,0.0775855,-8.62142e-05,4.21855e-08,-5.99056e-12,78708,30.4467], Tmin=(100,'K'), Tmax=(889.614,'K')), NASAPolynomial(coeffs=[18.4581,0.0119977,-3.04332e-06,4.05087e-10,-2.35403e-14,74825.7,-58.8619], Tmin=(889.614,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
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
    label = '[CH]=[C]C1C=C=COO1(25214)',
    structure = SMILES('[CH]=[C]C1C=C=COO1'),
    E0 = (701.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.909029,0.11813,-0.000180299,1.36136e-07,-3.84464e-11,84546.9,15.7001], Tmin=(100,'K'), Tmax=(662.628,'K')), NASAPolynomial(coeffs=[15.622,0.0321699,-1.7019e-05,3.36069e-09,-2.36373e-13,82052.5,-59.5405], Tmin=(662.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(701.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=COOC#C[CH2](25215)',
    structure = SMILES('[CH]=C=COOC#C[CH2]'),
    E0 = (703.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055,350,500,795,815,3000,3100,440,815,1455,1000,2100,2250,500,550],'cm^-1')),
        HinderedRotor(inertia=(1.53662,'amu*angstrom^2'), symmetry=1, barrier=(35.33,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53959,'amu*angstrom^2'), symmetry=1, barrier=(35.3981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53852,'amu*angstrom^2'), symmetry=1, barrier=(35.3736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53681,'amu*angstrom^2'), symmetry=1, barrier=(35.3343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580859,0.0741155,-8.47249e-05,4.92384e-08,-1.12825e-11,84789,29.8764], Tmin=(100,'K'), Tmax=(1066.41,'K')), NASAPolynomial(coeffs=[15.0292,0.0199201,-8.49261e-06,1.5807e-09,-1.09759e-13,81707.5,-40.7572], Tmin=(1066.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(703.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtOs) + radical(C=C=CJ) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]OOC=C=C(25216)',
    structure = SMILES('[CH]=C=[C]OOC=C=C'),
    E0 = (738.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,563.333,586.667,610,1970,2140,350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.17856,'amu*angstrom^2'), symmetry=1, barrier=(27.0975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17794,'amu*angstrom^2'), symmetry=1, barrier=(27.0831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17811,'amu*angstrom^2'), symmetry=1, barrier=(27.0871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615738,0.0763459,-9.79251e-05,6.5491e-08,-1.73232e-11,88944.1,31.6574], Tmin=(100,'K'), Tmax=(926.812,'K')), NASAPolynomial(coeffs=[13.3094,0.0215594,-9.2522e-06,1.7051e-09,-1.16785e-13,86591.3,-28.6171], Tmin=(926.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
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
    label = '[CH]=C=CO[O](20803)',
    structure = SMILES('[CH]=C=CO[O]'),
    E0 = (404.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,540,610,2055,180,180,180,2867.84],'cm^-1')),
        HinderedRotor(inertia=(0.0105329,'amu*angstrom^2'), symmetry=1, barrier=(9.09612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.1,'J/mol'), sigma=(5.7752,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.83 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13751,0.0383774,-4.88548e-05,3.09451e-08,-7.37348e-12,48769,18.1869], Tmin=(100,'K'), Tmax=(1167.31,'K')), NASAPolynomial(coeffs=[10.2545,0.00577971,-8.1986e-07,1.20449e-12,5.54427e-15,47199.8,-20.8326], Tmin=(1167.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(ROOJ)"""),
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
    label = '[CH]=C=[C]OOC=C=[CH](25217)',
    structure = SMILES('[CH]=C=[C]OOC=C=[CH]'),
    E0 = (893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,1685,370,540,563.333,586.667,610,1970,2140,350,500,795,815,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.37622,'amu*angstrom^2'), symmetry=1, barrier=(31.6421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37545,'amu*angstrom^2'), symmetry=1, barrier=(31.6244,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.606953,0.0774244,-0.00010507,6.81796e-08,-1.54003e-11,107523,32.3249], Tmin=(100,'K'), Tmax=(744.123,'K')), NASAPolynomial(coeffs=[14.1813,0.0163786,-6.04712e-06,9.95608e-10,-6.25939e-14,105173,-31.3695], Tmin=(744.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=COOC1[C]=C1(25218)',
    structure = SMILES('[CH]=C=COOC1[C]=C1'),
    E0 = (796.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.453136,0.0728961,-7.84561e-05,4.16505e-08,-8.617e-12,95975.9,29.3597], Tmin=(100,'K'), Tmax=(1184.57,'K')), NASAPolynomial(coeffs=[17.2517,0.0161712,-6.62589e-06,1.22475e-09,-8.52117e-14,91996.1,-54.5287], Tmin=(1184.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(796.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C1[CH]OOC=C=C1(25095)',
    structure = SMILES('[CH]=C1[CH]OOC=C=C1'),
    E0 = (624.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16928,0.0424995,2.41498e-05,-7.08933e-08,3.21368e-11,75249.3,19.2385], Tmin=(100,'K'), Tmax=(963.062,'K')), NASAPolynomial(coeffs=[20.2661,0.0117744,-3.67841e-06,7.61333e-10,-6.40543e-14,69317.6,-83.8732], Tmin=(963.062,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]COO[C]=C=[CH](25219)',
    structure = SMILES('[CH]=[C]COO[C]=C=[CH]'),
    E0 = (1049.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,2750,2850,1437.5,1250,1305,750,350,540,610,2055,350,500,795,815,1670,1700,300,440],'cm^-1')),
        HinderedRotor(inertia=(3.30403,'amu*angstrom^2'), symmetry=1, barrier=(75.9661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.496742,'amu*angstrom^2'), symmetry=1, barrier=(11.4211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.497074,'amu*angstrom^2'), symmetry=1, barrier=(11.4287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.47361,'amu*angstrom^2'), symmetry=1, barrier=(102.857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464519,0.0866127,-0.000143962,1.2689e-07,-4.31194e-11,126347,34.6309], Tmin=(100,'K'), Tmax=(863.512,'K')), NASAPolynomial(coeffs=[8.82933,0.0299037,-1.42527e-05,2.66196e-09,-1.79483e-13,125572,-0.619553], Tmin=(863.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1049.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=[C]OO[CH]C=[CH](25220)',
    structure = SMILES('[CH]=C=[C]OO[CH]C=[CH]'),
    E0 = (928.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3025,407.5,1350,352.5,540,610,2055,350,500,795,815,1685,370,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.805079,'amu*angstrom^2'), symmetry=1, barrier=(25.6273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18252,'amu*angstrom^2'), symmetry=1, barrier=(69.4735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18252,'amu*angstrom^2'), symmetry=1, barrier=(69.4735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805088,'amu*angstrom^2'), symmetry=1, barrier=(25.6273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753571,0.0703175,-7.65048e-05,4.26494e-08,-9.44534e-12,111848,33.7208], Tmin=(100,'K'), Tmax=(1097.36,'K')), NASAPolynomial(coeffs=[14.2762,0.0210256,-9.12651e-06,1.71556e-09,-1.1977e-13,108880,-32.7746], Tmin=(1097.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(928.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_P) + radical(C=CCJO) + radical(C=CJO)"""),
)

species(
    label = 'C1=COC=1(22275)',
    structure = SMILES('C1=COC=1'),
    E0 = (388.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17882,0.00883829,2.87655e-05,-4.70958e-08,1.95464e-11,46730.5,10.2865], Tmin=(100,'K'), Tmax=(941.322,'K')), NASAPolynomial(coeffs=[9.87744,0.00380291,-5.45265e-07,1.04211e-10,-1.15627e-14,44431.4,-27.139], Tmin=(941.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[CH]=C=COC(=[CH])C=O(22650)',
    structure = SMILES('[CH]=C=COC(=[CH])C=O'),
    E0 = (390.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,2782.5,750,1395,475,1775,1000,540,610,2055,350,440,435,1725,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.06795,'amu*angstrom^2'), symmetry=1, barrier=(24.5544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06925,'amu*angstrom^2'), symmetry=1, barrier=(24.5842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06896,'amu*angstrom^2'), symmetry=1, barrier=(24.5775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3842.2,'J/mol'), sigma=(6.0658,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=600.14 K, Pc=39.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475741,0.0784617,-9.66611e-05,5.96333e-08,-1.44515e-11,47118.1,27.4708], Tmin=(100,'K'), Tmax=(1010.53,'K')), NASAPolynomial(coeffs=[15.4177,0.0193148,-8.86265e-06,1.70915e-09,-1.2092e-13,44098.3,-44.7715], Tmin=(1010.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C1=COO[CH]C1=[CH](25075)',
    structure = SMILES('[CH]C1=COO[CH]C1=[CH]'),
    E0 = (736.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15998,0.0354422,6.99442e-05,-1.22296e-07,5.00123e-11,88664.2,20.8871], Tmin=(100,'K'), Tmax=(971.426,'K')), NASAPolynomial(coeffs=[21.652,0.018112,-6.82771e-06,1.44184e-09,-1.17518e-13,81519.3,-93.6644], Tmin=(971.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(736.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=CCJO) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]#C[CH]OOC=C=[CH](25221)',
    structure = SMILES('[C]#C[CH]OOC=C=[CH]'),
    E0 = (1014.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,540,610,2055,350,500,795,815,2175,525,180],'cm^-1')),
        HinderedRotor(inertia=(2.48929,'amu*angstrom^2'), symmetry=1, barrier=(57.2336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.49085,'amu*angstrom^2'), symmetry=1, barrier=(57.2696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.49114,'amu*angstrom^2'), symmetry=1, barrier=(57.2761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.49131,'amu*angstrom^2'), symmetry=1, barrier=(57.2801,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.308248,0.0856905,-0.000133651,1.05163e-07,-3.14644e-11,122130,29.8399], Tmin=(100,'K'), Tmax=(948.512,'K')), NASAPolynomial(coeffs=[12.8563,0.019209,-7.06392e-06,1.1125e-09,-6.57535e-14,120360,-26.817], Tmin=(948.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1014.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOOC) + radical(Acetyl) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]C1OOC1C#C(25222)',
    structure = SMILES('[CH]=[C]C1OOC1C#C'),
    E0 = (775.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.592895,0.0718294,-7.3153e-05,3.62151e-08,-6.48947e-12,93476.6,32.0817], Tmin=(100,'K'), Tmax=(1691.39,'K')), NASAPolynomial(coeffs=[17.2559,0.00850291,1.73376e-06,-6.82617e-10,5.52293e-14,90459.1,-54.4805], Tmin=(1691.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(775.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxetane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[C]#CCOOC=C=[CH](25223)',
    structure = SMILES('[C]#CCOOC=C=[CH]'),
    E0 = (827.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,2175,525,180],'cm^-1')),
        HinderedRotor(inertia=(2.40024,'amu*angstrom^2'), symmetry=1, barrier=(55.1864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.39485,'amu*angstrom^2'), symmetry=1, barrier=(55.0623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.39605,'amu*angstrom^2'), symmetry=1, barrier=(55.0899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993438,'amu*angstrom^2'), symmetry=1, barrier=(22.8411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582351,0.0798485,-0.000101602,5.70354e-08,-7.01136e-12,99696.6,28.6001], Tmin=(100,'K'), Tmax=(676.547,'K')), NASAPolynomial(coeffs=[13.5258,0.0203042,-7.23514e-06,1.14827e-09,-6.96519e-14,97556.6,-31.659], Tmin=(676.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(827.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(Acetyl) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=[C]OOCC#C(25224)',
    structure = SMILES('[CH]=C=[C]OOCC#C'),
    E0 = (730.536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,540,610,2055,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,1685,370,2175,525],'cm^-1')),
        HinderedRotor(inertia=(0.729227,'amu*angstrom^2'), symmetry=1, barrier=(16.7664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43947,'amu*angstrom^2'), symmetry=1, barrier=(56.0883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43961,'amu*angstrom^2'), symmetry=1, barrier=(56.0915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43906,'amu*angstrom^2'), symmetry=1, barrier=(56.0788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510125,0.0817759,-0.000122761,9.84033e-08,-3.07166e-11,87984.1,31.2081], Tmin=(100,'K'), Tmax=(893.95,'K')), NASAPolynomial(coeffs=[10.8707,0.0251489,-1.05143e-05,1.84521e-09,-1.19993e-13,86542.1,-15.3193], Tmin=(893.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(730.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[C]#C[CH]OOC=C=C(25225)',
    structure = SMILES('[C]#C[CH]OOC=C=C'),
    E0 = (859.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,350,500,795,815,540,610,2055,2950,3100,1380,975,1025,1650,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.87176,'amu*angstrom^2'), symmetry=1, barrier=(66.0273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87895,'amu*angstrom^2'), symmetry=1, barrier=(66.1928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36045,'amu*angstrom^2'), symmetry=1, barrier=(31.2794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.86812,'amu*angstrom^2'), symmetry=1, barrier=(65.9437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375613,0.0837213,-0.000122244,9.47836e-08,-2.88379e-11,103549,28.9737], Tmin=(100,'K'), Tmax=(878.238,'K')), NASAPolynomial(coeffs=[12.1057,0.0241978,-1.01645e-05,1.79842e-09,-1.18061e-13,101724,-24.7552], Tmin=(878.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(859.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOOC) + radical(Acetyl)"""),
)

species(
    label = '[CH]C1=COOC1C#C(25226)',
    structure = SMILES('[CH]C1=COOC1C#C'),
    E0 = (561.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.91927,0.0531127,-4.61569e-06,-3.93229e-08,2.09071e-11,67616.1,26.3308], Tmin=(100,'K'), Tmax=(951.73,'K')), NASAPolynomial(coeffs=[17.784,0.0180396,-5.773e-06,1.01949e-09,-7.41683e-14,62784.3,-62.7174], Tmin=(951.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxolene) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]1=COOC=[C]C=C1(25144)',
    structure = SMILES('[C]1=COOC=[C]C=C1'),
    E0 = (569.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72199,0.0246986,7.91798e-05,-1.34129e-07,5.73625e-11,68558,18.041], Tmin=(100,'K'), Tmax=(920.683,'K')), NASAPolynomial(coeffs=[20.9522,0.00678973,1.41697e-06,-3.85229e-10,1.94154e-14,62235.1,-88.2525], Tmin=(920.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(C=CJC=C) + radical(C=CJC=C)"""),
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
    label = '[CH]OOC=C=[CH](23157)',
    structure = SMILES('[CH]OOC=C=[CH]'),
    E0 = (713.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,500,795,815,3010,987.5,1337.5,450,1655,311.097,851.355],'cm^-1')),
        HinderedRotor(inertia=(0.594343,'amu*angstrom^2'), symmetry=1, barrier=(13.6651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59463,'amu*angstrom^2'), symmetry=1, barrier=(13.6717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.84824,'amu*angstrom^2'), symmetry=1, barrier=(88.4785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31984,0.0585742,-7.51673e-05,4.76862e-08,-1.1737e-11,85949.5,22.4565], Tmin=(100,'K'), Tmax=(1002.18,'K')), NASAPolynomial(coeffs=[13.0987,0.0115612,-4.80092e-06,8.7714e-10,-6.01566e-14,83588.6,-34.3952], Tmin=(1002.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CH2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#C[CH]OOCC#C(25227)',
    structure = SMILES('[C]#C[CH]OOCC#C'),
    E0 = (851.913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,350,500,795,815,2100,2250,500,550,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.76967,'amu*angstrom^2'), symmetry=1, barrier=(63.6801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.77021,'amu*angstrom^2'), symmetry=1, barrier=(63.6926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.77054,'amu*angstrom^2'), symmetry=1, barrier=(63.7002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0301167,'amu*angstrom^2'), symmetry=1, barrier=(63.7694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.77202,'amu*angstrom^2'), symmetry=1, barrier=(63.7341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375182,0.087837,-0.000142122,1.20707e-07,-3.89818e-11,102584,28.1518], Tmin=(100,'K'), Tmax=(918.278,'K')), NASAPolynomial(coeffs=[9.40782,0.028241,-1.16938e-05,2.00258e-09,-1.26641e-13,101779,-10.0066], Tmin=(918.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(851.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CtOsHH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl) + radical(CCsJOOC)"""),
)

species(
    label = 'C#CC1OOC1C#C(22654)',
    structure = SMILES('C#CC1OOC1C#C'),
    E0 = (456.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14229,0.0737762,-7.48832e-05,3.68568e-08,-6.45608e-12,55140,30.1195], Tmin=(100,'K'), Tmax=(1781.49,'K')), NASAPolynomial(coeffs=[15.4046,0.00897836,2.95279e-06,-9.8142e-10,7.61921e-14,53631.3,-46.9521], Tmin=(1781.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + ring(12dioxetane)"""),
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
    E0 = (653.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (778.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (870.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (983.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (653.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (912.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1104.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (803.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (756.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1088.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (968.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (743.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (967.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (736.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1226.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (775.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (653.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (980.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (774.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1042.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (809.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (730.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1300.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (935.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (661.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['C#CC=O(21959)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['[CH]=[C]C1C=C=COO1(25214)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2e+11,'s^-1'), n=0.21, Ea=(125.52,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R7_MMSR;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C=COOC#C[CH2](25215)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.65028e+09,'s^-1'), n=1.32317, Ea=(166.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] + [R3H;Cd_rad_out_singleNd;XH_out] for rate rule [R3H;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C=[C]OOC=C=C(25216)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_double;Cd_H_out_singleH] for rate rule [R6H;Cd_rad_out_double;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C=C[O](8556)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(905000,'m^3/(mol*s)'), n=0, Ea=(113.507,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;O_sec_rad] for rate rule [O_rad/OneDe;O_rad/OneDe]
Euclidian distance = 1.41421356237
family: R_Recombination
Ea raised from 0.0 to 113.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C=CO[O](20803)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.15767e+07,'m^3/(mol*s)'), n=0.0716491, Ea=(15.4197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad/NonDe;Y_rad] for rate rule [O_rad/NonDe;Cd_allenic]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(8)', '[CH]=C=[C]OOC=C=[CH](25217)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['[CH]=C=COOC1[C]=C1(25218)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.4708e+12,'s^-1'), n=-0.1205, Ea=(149.998,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['[CH]=C1[CH]OOC=C=C1(25095)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.70414e+11,'s^-1'), n=0.0116667, Ea=(102.85,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra;radadd_intra_cdsingleH] + [R7_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]COO[C]=C=[CH](25219)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C=[C]OO[CH]C=[CH](25220)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['[CH]=C=C[O](8556)', 'C1=COC=1(22275)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.22709e+11,'s^-1'), n=0, Ea=(89.9403,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OOR] for rate rule [R3OO;Cd_pri_rad_in;OOR]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['[CH]=C=COC(=[CH])C=O(22650)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(14080,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C1=COO[CH]C1=[CH](25075)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[C]#C[CH]OOC=C=[CH](25221)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['[CH]=[C]C1OOC1C#C(25222)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(413769,'s^-1'), n=1.87624, Ea=(122.367,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 119.4 to 122.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C=C[O](8556)', 'C#CC=O(21959)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(372.961,'m^3/(mol*s)'), n=1.215, Ea=(299.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Od_CO-CtH;YJ] for rate rule [Od_CO-CtH;O_rad/OneDe]
Euclidian distance = 4.0
family: R_Addition_MultipleBond
Ea raised from 295.2 to 299.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]#CCOOC=C=[CH](25223)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.9263e+09,'s^-1'), n=1.08337, Ea=(153.033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C=[C]OOCC#C(25224)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_H/Ct]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[C]#C[CH]OOC=C=C(25225)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(640643,'s^-1'), n=2.07799, Ea=(182.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Y_rad_out;Cd_H_out_singleH] for rate rule [R8Hall;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['[CH]C1=COOC1C#C(25226)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.3687e+16,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['[C]1=COOC=[C]C=C1(25144)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.53286e+10,'s^-1'), n=0.161, Ea=(76.7659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;triplebond_intra_H;radadd_intra_cdsingleH] for rate rule [R8_linear;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]#C(5143)', '[CH]OOC=C=[CH](23157)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[C]#C[CH]OOCC#C(25227)'],
    products = ['[CH]=C=COOC=C=[CH](22648)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(49528.1,'s^-1'), n=1.95205, Ea=(83.4098,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C=COOC=C=[CH](22648)'],
    products = ['C#CC1OOC1C#C(22654)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

network(
    label = '4951',
    isomers = [
        '[CH]=C=COOC=C=[CH](22648)',
    ],
    reactants = [
        ('C#CC=O(21959)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4951',
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

