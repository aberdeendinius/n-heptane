species(
    label = '[CH2]C([CH]O)=C[C]=C(20132)',
    structure = SMILES('[CH2]C([CH]O)=C[C]=C'),
    E0 = (300.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,473.899,475.038],'cm^-1')),
        HinderedRotor(inertia=(0.177816,'amu*angstrom^2'), symmetry=1, barrier=(28.5272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.437639,'amu*angstrom^2'), symmetry=1, barrier=(69.6151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434448,'amu*angstrom^2'), symmetry=1, barrier=(69.5485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434621,'amu*angstrom^2'), symmetry=1, barrier=(69.5933,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970545,0.0546111,-1.36565e-05,-2.84809e-08,1.73025e-11,36207.9,26.9384], Tmin=(100,'K'), Tmax=(930.571,'K')), NASAPolynomial(coeffs=[16.3852,0.0189076,-5.35797e-06,8.58868e-10,-5.90609e-14,32016.1,-53.4277], Tmin=(930.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CCJO)"""),
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
    label = '[CH2]C(C=C=C)=C[O](20110)',
    structure = SMILES('[CH2]C(C=C=C)=C[O]'),
    E0 = (281.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43665,'amu*angstrom^2'), symmetry=1, barrier=(33.0315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.437,'amu*angstrom^2'), symmetry=1, barrier=(33.0395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612303,0.058359,-1.49777e-05,-3.7588e-08,2.29723e-11,33942.8,23.1652], Tmin=(100,'K'), Tmax=(929.87,'K')), NASAPolynomial(coeffs=[21.8819,0.00821482,-7.92924e-07,6.55475e-11,-8.55754e-15,28199.5,-87.5145], Tmin=(929.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C#CC([CH2])=CO(28959)',
    structure = SMILES('[CH2]C#CC([CH2])=CO'),
    E0 = (287.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2100,2250,500,550,350,440,435,1725,219.105],'cm^-1')),
        HinderedRotor(inertia=(0.795362,'amu*angstrom^2'), symmetry=1, barrier=(27.382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.797018,'amu*angstrom^2'), symmetry=1, barrier=(27.4175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800078,'amu*angstrom^2'), symmetry=1, barrier=(27.3794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800028,'amu*angstrom^2'), symmetry=1, barrier=(27.3622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759286,0.056069,-1.39305e-05,-3.71547e-08,2.28159e-11,34671,24.712], Tmin=(100,'K'), Tmax=(920.536,'K')), NASAPolynomial(coeffs=[20.9562,0.00804506,-4.28103e-07,-3.90985e-11,4.41362e-16,29269,-80.1992], Tmin=(920.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Allyl_P) + radical(Propargyl)"""),
)

species(
    label = 'C#CC=C([CH2])[CH]O(27790)',
    structure = SMILES('C#CC=C([CH2])[CH]O'),
    E0 = (281.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2175,525,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,750,770,3400,2100,504.061],'cm^-1')),
        HinderedRotor(inertia=(0.140798,'amu*angstrom^2'), symmetry=1, barrier=(25.3859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.412031,'amu*angstrom^2'), symmetry=1, barrier=(74.2886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140799,'amu*angstrom^2'), symmetry=1, barrier=(25.3859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.412031,'amu*angstrom^2'), symmetry=1, barrier=(74.2887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02852,0.0549948,-2.64598e-05,-1.09554e-08,9.88006e-12,33937.2,24.1194], Tmin=(100,'K'), Tmax=(961.409,'K')), NASAPolynomial(coeffs=[15.8959,0.0174461,-5.80127e-06,1.01785e-09,-7.19314e-14,29955.1,-52.8643], Tmin=(961.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(C=CCJO)"""),
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
    label = 'C=[C]C=C=CO(28960)',
    structure = SMILES('C=[C]C=C=CO'),
    E0 = (224.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.43509,'amu*angstrom^2'), symmetry=1, barrier=(32.9956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43908,'amu*angstrom^2'), symmetry=1, barrier=(33.0874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.220575,0.0681799,-7.31183e-05,3.65318e-08,-6.60853e-12,27198.7,23.4728], Tmin=(100,'K'), Tmax=(1643.39,'K')), NASAPolynomial(coeffs=[19.2481,0.00307808,2.47256e-06,-6.92169e-10,5.19901e-14,23192,-72.8455], Tmin=(1643.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]O(5471)',
    structure = SMILES('[CH]O'),
    E0 = (205.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,402.686,3356.18],'cm^-1')),
        HinderedRotor(inertia=(0.0105042,'amu*angstrom^2'), symmetry=1, barrier=(23.1306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76003,0.0029575,8.86344e-06,-1.3392e-08,5.33433e-12,24775.7,6.76105], Tmin=(100,'K'), Tmax=(943.117,'K')), NASAPolynomial(coeffs=[5.07489,0.00326005,-9.68482e-07,1.67779e-10,-1.21779e-14,24266.2,-0.891576], Tmin=(943.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = 'C=[C]C=C=C(19283)',
    structure = SMILES('C=[C]C=C=C'),
    E0 = (433.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.5987,'amu*angstrom^2'), symmetry=1, barrier=(36.7571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98554,0.0388285,-2.1811e-05,-2.28683e-09,4.80018e-12,52215,14.9887], Tmin=(100,'K'), Tmax=(940.235,'K')), NASAPolynomial(coeffs=[10.8095,0.0147077,-4.73738e-06,7.85957e-10,-5.27182e-14,49962.5,-30.1921], Tmin=(940.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(=C[C]=C)C[O](20752)',
    structure = SMILES('[CH2]C(=C[C]=C)C[O]'),
    E0 = (408.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180,971.914,2516.91],'cm^-1')),
        HinderedRotor(inertia=(3.02026,'amu*angstrom^2'), symmetry=1, barrier=(69.4418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.02122,'amu*angstrom^2'), symmetry=1, barrier=(69.4638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.02036,'amu*angstrom^2'), symmetry=1, barrier=(69.4441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39332,0.0596304,-5.38485e-05,3.06909e-08,-7.70845e-12,49218.4,25.6931], Tmin=(100,'K'), Tmax=(929.754,'K')), NASAPolynomial(coeffs=[7.01462,0.035446,-1.48304e-05,2.71308e-09,-1.85433e-13,48173.1,-1.01694], Tmin=(929.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CCOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=[C]C([CH2])=CO(28961)',
    structure = SMILES('[CH2]C=[C]C([CH2])=CO'),
    E0 = (280.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.362915,0.063283,-1.47261e-05,-4.64416e-08,2.90811e-11,33830.5,24.7542], Tmin=(100,'K'), Tmax=(894.878,'K')), NASAPolynomial(coeffs=[23.3503,0.00703103,1.62215e-06,-5.55649e-10,4.04941e-14,27854.5,-93.9962], Tmin=(894.878,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C=CC([CH2])=CO(28744)',
    structure = SMILES('[CH]C=CC([CH2])=CO'),
    E0 = (333.696,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.98561,'amu*angstrom^2'), symmetry=1, barrier=(45.653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98367,'amu*angstrom^2'), symmetry=1, barrier=(45.6085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98622,'amu*angstrom^2'), symmetry=1, barrier=(45.667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98529,'amu*angstrom^2'), symmetry=1, barrier=(45.6457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.161647,0.0655937,-8.8067e-06,-5.31581e-08,3.05363e-11,40290.1,25.4285], Tmin=(100,'K'), Tmax=(916.865,'K')), NASAPolynomial(coeffs=[23.9154,0.0114095,-1.05539e-06,2.55007e-11,-3.46247e-15,33856,-98.4417], Tmin=(916.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=[C][C]=C)CO(28962)',
    structure = SMILES('[CH2]C(=[C][C]=C)CO'),
    E0 = (381.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,180,1456.93],'cm^-1')),
        HinderedRotor(inertia=(3.81918,'amu*angstrom^2'), symmetry=1, barrier=(87.8104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.8203,'amu*angstrom^2'), symmetry=1, barrier=(87.8363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.82097,'amu*angstrom^2'), symmetry=1, barrier=(87.8517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.82352,'amu*angstrom^2'), symmetry=1, barrier=(87.9103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940649,0.0689112,-6.96834e-05,3.54076e-08,-4.96313e-12,46023,26.0626], Tmin=(100,'K'), Tmax=(749.019,'K')), NASAPolynomial(coeffs=[10.513,0.0285328,-1.0331e-05,1.72593e-09,-1.11245e-13,44287.7,-19.3633], Tmin=(749.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=[C]C(C)=CO(28963)',
    structure = SMILES('[CH2][C]=[C]C(C)=CO'),
    E0 = (366.405,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.33202,'amu*angstrom^2'), symmetry=1, barrier=(30.6257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33166,'amu*angstrom^2'), symmetry=1, barrier=(30.6174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33416,'amu*angstrom^2'), symmetry=1, barrier=(30.6749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34329,'amu*angstrom^2'), symmetry=1, barrier=(30.8849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.504436,0.080264,-8.30004e-05,4.23014e-08,-8.05267e-12,44246.7,28.0576], Tmin=(100,'K'), Tmax=(1495.61,'K')), NASAPolynomial(coeffs=[20.2675,0.0116046,-9.96001e-07,-1.10668e-10,1.60489e-14,39499.1,-75.6161], Tmin=(1495.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CC(C)=C[O](20112)',
    structure = SMILES('[CH2][C]=CC(C)=C[O]'),
    E0 = (308.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23315,'amu*angstrom^2'), symmetry=1, barrier=(28.3525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23613,'amu*angstrom^2'), symmetry=1, barrier=(28.421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2353,'amu*angstrom^2'), symmetry=1, barrier=(28.4021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.81108,0.0570286,-1.65239e-05,-2.85268e-08,1.79064e-11,37275.6,25.1487], Tmin=(100,'K'), Tmax=(934.525,'K')), NASAPolynomial(coeffs=[18.0418,0.0163698,-4.38014e-06,7.02768e-10,-4.98156e-14,32610.1,-64.5442], Tmin=(934.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC([CH2])=C[O](20116)',
    structure = SMILES('[CH2]C=CC([CH2])=C[O]'),
    E0 = (222.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.846573,0.0498828,1.7124e-05,-7.1807e-08,3.53018e-11,26895.7,24.829], Tmin=(100,'K'), Tmax=(924.874,'K')), NASAPolynomial(coeffs=[21.4792,0.0108938,-1.13186e-06,9.16153e-11,-1.05796e-14,20930.3,-84.7177], Tmin=(924.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]C=C([CH2])CO(28756)',
    structure = SMILES('[CH]=[C]C=C([CH2])CO'),
    E0 = (429.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,251.972],'cm^-1')),
        HinderedRotor(inertia=(2.28833,'amu*angstrom^2'), symmetry=1, barrier=(103.183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28724,'amu*angstrom^2'), symmetry=1, barrier=(103.167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136005,'amu*angstrom^2'), symmetry=1, barrier=(6.12406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28774,'amu*angstrom^2'), symmetry=1, barrier=(103.172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716025,0.068894,-6.93782e-05,3.84628e-08,-8.56651e-12,51820.4,27.6562], Tmin=(100,'K'), Tmax=(1091.79,'K')), NASAPolynomial(coeffs=[12.7521,0.0247965,-8.79214e-06,1.46723e-09,-9.50619e-14,49192.2,-31.4678], Tmin=(1091.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#C[CH][C](C)[CH]O(28745)',
    structure = SMILES('C#C[CH][C](C)[CH]O'),
    E0 = (412.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,750,770,3400,2100,3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,3615,1277.5,1000,299.775,2872.1],'cm^-1')),
        HinderedRotor(inertia=(0.280516,'amu*angstrom^2'), symmetry=1, barrier=(17.8884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18655,'amu*angstrom^2'), symmetry=1, barrier=(75.6658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18655,'amu*angstrom^2'), symmetry=1, barrier=(75.6658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18655,'amu*angstrom^2'), symmetry=1, barrier=(75.6658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0129262,'amu*angstrom^2'), symmetry=1, barrier=(75.6659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.780936,0.0714448,-7.63767e-05,4.11304e-08,-7.07113e-12,49695.9,25.9824], Tmin=(100,'K'), Tmax=(783.615,'K')), NASAPolynomial(coeffs=[11.9881,0.0251803,-8.76432e-06,1.42987e-09,-9.086e-14,47603.4,-27.497], Tmin=(783.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJ(C)CO) + radical(Sec_Propargyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][C]=CC([CH2])=C[O](20121)',
    structure = SMILES('[CH2][C]=CC([CH2])=C[O]'),
    E0 = (460.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.67908,'amu*angstrom^2'), symmetry=1, barrier=(38.6054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67688,'amu*angstrom^2'), symmetry=1, barrier=(38.5547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67364,'amu*angstrom^2'), symmetry=1, barrier=(38.4803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.81984,0.0544109,-7.43188e-06,-4.27957e-08,2.44201e-11,55498.8,25.3986], Tmin=(100,'K'), Tmax=(923.115,'K')), NASAPolynomial(coeffs=[20.2652,0.0104351,-1.43316e-06,1.46181e-10,-1.23669e-14,50192.4,-76.1553], Tmin=(923.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CC=CCJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=[C][C]=C)[CH]O(28964)',
    structure = SMILES('[CH2]C(=[C][C]=C)[CH]O'),
    E0 = (499.047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,417.894,417.905],'cm^-1')),
        HinderedRotor(inertia=(0.447279,'amu*angstrom^2'), symmetry=1, barrier=(55.4633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447532,'amu*angstrom^2'), symmetry=1, barrier=(55.4623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447438,'amu*angstrom^2'), symmetry=1, barrier=(55.4623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447622,'amu*angstrom^2'), symmetry=1, barrier=(55.4634,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381086,0.0671605,-6.36613e-05,3.14372e-08,-5.99271e-12,60162.3,28.6519], Tmin=(100,'K'), Tmax=(1428.54,'K')), NASAPolynomial(coeffs=[15.951,0.0170751,-4.25709e-06,5.35022e-10,-2.82821e-14,56376,-49.6995], Tmin=(1428.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CCJO) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C]C=C([CH2])[CH]O(21903)',
    structure = SMILES('[CH]=[C]C=C([CH2])[CH]O'),
    E0 = (547.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,433.307],'cm^-1')),
        HinderedRotor(inertia=(0.50806,'amu*angstrom^2'), symmetry=1, barrier=(67.0981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507427,'amu*angstrom^2'), symmetry=1, barrier=(67.0725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224943,'amu*angstrom^2'), symmetry=1, barrier=(29.9272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.501106,'amu*angstrom^2'), symmetry=1, barrier=(67.0762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902672,0.0580589,-3.00243e-05,-1.15182e-08,1.15877e-11,65927.4,27.5855], Tmin=(100,'K'), Tmax=(922.595,'K')), NASAPolynomial(coeffs=[16.6204,0.0160842,-4.32983e-06,6.61554e-10,-4.43179e-14,61913.4,-53.0132], Tmin=(922.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]1C(O)C1[C]=C(28965)',
    structure = SMILES('[CH2][C]1C(O)C1[C]=C'),
    E0 = (504.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47421,0.0413066,1.59126e-05,-5.50344e-08,2.5742e-11,60784.7,26.7766], Tmin=(100,'K'), Tmax=(937.62,'K')), NASAPolynomial(coeffs=[15.2232,0.0185855,-5.22547e-06,8.69742e-10,-6.2386e-14,56626.9,-47.092], Tmin=(937.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJ(C)CO) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C]C1C[C]1[CH]O(28966)',
    structure = SMILES('C=[C]C1C[C]1[CH]O'),
    E0 = (492.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20795,0.0503224,-1.29909e-05,-2.14681e-08,1.248e-11,59330.9,26.5379], Tmin=(100,'K'), Tmax=(980.657,'K')), NASAPolynomial(coeffs=[14.5517,0.0213645,-7.6556e-06,1.38946e-09,-9.88223e-14,55489,-43.8217], Tmin=(980.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(492.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCsJOH) + radical(Cds_S) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C1([CH]O)[CH]C1=C(28801)',
    structure = SMILES('[CH2]C1([CH]O)[CH]C1=C'),
    E0 = (489.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.742912,0.0613986,-3.82196e-05,1.71807e-09,4.64236e-12,58947.4,24.2138], Tmin=(100,'K'), Tmax=(1019.25,'K')), NASAPolynomial(coeffs=[16.3861,0.0203806,-7.83742e-06,1.45696e-09,-1.03827e-13,54700.3,-56.7453], Tmin=(1019.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(CCsJOH) + radical(Allyl_S) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C#CC(=C)CO(28967)',
    structure = SMILES('[CH2]C#CC(=C)CO'),
    E0 = (191.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38118,0.0557923,-4.32113e-05,1.81607e-08,-3.19526e-12,23108.7,26.0934], Tmin=(100,'K'), Tmax=(1315.54,'K')), NASAPolynomial(coeffs=[10.5217,0.0279997,-1.15217e-05,2.10145e-09,-1.43416e-13,20703.8,-20.5111], Tmin=(1315.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl)"""),
)

species(
    label = '[CH2]C#CC(C)=CO(28968)',
    structure = SMILES('[CH2]C#CC(C)=CO'),
    E0 = (135.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751093,0.0586806,-2.30048e-05,-2.29017e-08,1.63048e-11,16447.8,24.4601], Tmin=(100,'K'), Tmax=(931.448,'K')), NASAPolynomial(coeffs=[18.7274,0.0139889,-3.38028e-06,5.18711e-10,-3.71084e-14,11689,-68.5579], Tmin=(931.448,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl)"""),
)

species(
    label = 'C=[C]C=C(C)C=O(20136)',
    structure = SMILES('C=[C]C=C(C)C=O'),
    E0 = (126.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32095,0.0622151,-5.43578e-05,2.63e-08,-5.44411e-12,15258.2,21.6455], Tmin=(100,'K'), Tmax=(1114.4,'K')), NASAPolynomial(coeffs=[9.39329,0.0332398,-1.53557e-05,2.96725e-09,-2.09616e-13,13459.1,-18.173], Tmin=(1114.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CC(=C)C=O(20126)',
    structure = SMILES('[CH2]C=CC(=C)C=O'),
    E0 = (61.9979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.999258,0.0532226,-1.35446e-05,-2.31093e-08,1.30285e-11,7575.91,23.3687], Tmin=(100,'K'), Tmax=(1005.29,'K')), NASAPolynomial(coeffs=[16.1789,0.0212779,-8.33634e-06,1.59218e-09,-1.16128e-13,3086.12,-57.0954], Tmin=(1005.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.9979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C=C[C]1CC1O(28969)',
    structure = SMILES('C=C=C[C]1CC1O'),
    E0 = (235.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53282,0.0377721,2.81109e-05,-6.66344e-08,2.89557e-11,28368.3,23.8989], Tmin=(100,'K'), Tmax=(961.876,'K')), NASAPolynomial(coeffs=[15.7687,0.0192193,-6.34502e-06,1.18022e-09,-8.8591e-14,23749.3,-54.002], Tmin=(961.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C=C1[CH]C(=C)C1O(28783)',
    structure = SMILES('C=C1[CH]C(=C)C1O'),
    E0 = (120.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96602,0.0188563,9.40572e-05,-1.42263e-07,5.758e-11,14533,23.071], Tmin=(100,'K'), Tmax=(941.501,'K')), NASAPolynomial(coeffs=[18.6659,0.0137364,-2.6666e-06,4.91115e-10,-4.57341e-14,8470.7,-71.9845], Tmin=(941.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C1[CH]C(=CO)C1(28772)',
    structure = SMILES('C=C1[CH]C(=CO)C1'),
    E0 = (89.0684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60896,0.0178283,0.00012341,-1.93679e-07,8.17328e-11,10830.9,21.829], Tmin=(100,'K'), Tmax=(917.099,'K')), NASAPolynomial(coeffs=[26.4967,-0.000599661,6.14806e-06,-1.2868e-09,7.76797e-14,2476.06,-116.749], Tmin=(917.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.0684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C)"""),
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
    label = '[CH]C([CH2])=C[C]=C(18814)',
    structure = SMILES('[CH]C([CH2])=C[C]=C'),
    E0 = (706.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,1685,370,3010,987.5,1337.5,450,1655,213.788,218.272,220.282,221.947],'cm^-1')),
        HinderedRotor(inertia=(1.47498,'amu*angstrom^2'), symmetry=1, barrier=(50.7819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49216,'amu*angstrom^2'), symmetry=1, barrier=(50.8696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53172,'amu*angstrom^2'), symmetry=1, barrier=(50.8545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3312.87,'J/mol'), sigma=(5.74688,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.46 K, Pc=39.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39389,0.0500704,-1.91482e-05,-1.09538e-08,8.38261e-12,85081.9,22.1753], Tmin=(100,'K'), Tmax=(944.304,'K')), NASAPolynomial(coeffs=[11.1102,0.0263476,-9.15979e-06,1.54654e-09,-1.03094e-13,82469.5,-28.2594], Tmin=(944.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=C[C]=C(19278)',
    structure = SMILES('[CH2][C]=C[C]=C'),
    E0 = (612.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(2.16813,'amu*angstrom^2'), symmetry=1, barrier=(49.8496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16803,'amu*angstrom^2'), symmetry=1, barrier=(49.8473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20158,0.0347752,-1.38743e-05,-8.03412e-09,6.49282e-12,73770.6,17.192], Tmin=(100,'K'), Tmax=(915.382,'K')), NASAPolynomial(coeffs=[9.16592,0.0169748,-5.40502e-06,8.73138e-10,-5.70753e-14,71966.4,-18.682], Tmin=(915.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=C[C]=CO(28970)',
    structure = SMILES('[CH2][C]=C[C]=CO'),
    E0 = (403.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1670,1700,300,440,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.59813,'amu*angstrom^2'), symmetry=1, barrier=(36.7441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59992,'amu*angstrom^2'), symmetry=1, barrier=(36.7852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59577,'amu*angstrom^2'), symmetry=1, barrier=(36.6898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0429243,0.0645806,-6.68019e-05,3.29802e-08,-5.91337e-12,48756.1,25.814], Tmin=(100,'K'), Tmax=(1667.6,'K')), NASAPolynomial(coeffs=[17.2779,0.00579288,1.58564e-06,-5.59288e-10,4.41954e-14,45376.5,-59.419], Tmin=(1667.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
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
    label = '[CH]C([CH2])=CO(19098)',
    structure = SMILES('[CH]C([CH2])=CO'),
    E0 = (280.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07603,'amu*angstrom^2'), symmetry=1, barrier=(47.7321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07698,'amu*angstrom^2'), symmetry=1, barrier=(47.7538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07646,'amu*angstrom^2'), symmetry=1, barrier=(47.7419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45652,0.0411178,9.02269e-06,-5.40178e-08,2.77613e-11,33836.9,18.3017], Tmin=(100,'K'), Tmax=(910.922,'K')), NASAPolynomial(coeffs=[18.0563,0.00742605,-4.82671e-08,-1.37232e-10,8.52061e-15,29186.3,-69.1609], Tmin=(910.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CC([CH2])=[C]O(28971)',
    structure = SMILES('[CH2][C]=CC([CH2])=[C]O'),
    E0 = (558.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,226.963],'cm^-1')),
        HinderedRotor(inertia=(0.709261,'amu*angstrom^2'), symmetry=1, barrier=(25.967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.705851,'amu*angstrom^2'), symmetry=1, barrier=(25.977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696119,'amu*angstrom^2'), symmetry=1, barrier=(26.0068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64168,'amu*angstrom^2'), symmetry=1, barrier=(61.1827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.289585,0.0736764,-7.3901e-05,3.59471e-08,-6.50941e-12,67362.7,31.11], Tmin=(100,'K'), Tmax=(1579.12,'K')), NASAPolynomial(coeffs=[19.9684,0.00951018,-7.41996e-07,-9.26983e-11,1.21963e-14,62567.1,-70.8046], Tmin=(1579.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C(C=[C][CH2])=CO(24296)',
    structure = SMILES('[CH]C(C=[C][CH2])=CO'),
    E0 = (538.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07941,'amu*angstrom^2'), symmetry=1, barrier=(47.8096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07942,'amu*angstrom^2'), symmetry=1, barrier=(47.8099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0794,'amu*angstrom^2'), symmetry=1, barrier=(47.8095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07939,'amu*angstrom^2'), symmetry=1, barrier=(47.8092,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416782,0.0636857,-1.85098e-05,-3.70321e-08,2.38442e-11,64861,26.2645], Tmin=(100,'K'), Tmax=(907.653,'K')), NASAPolynomial(coeffs=[21.1917,0.0128555,-1.80835e-06,1.33117e-10,-7.89505e-15,59412.2,-81.1903], Tmin=(907.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1([CH]O)C=[C]C1(28743)',
    structure = SMILES('[CH2]C1([CH]O)C=[C]C1'),
    E0 = (557.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.830846,0.0606483,-4.22893e-05,9.62022e-09,1.01413e-12,67229.4,26.6044], Tmin=(100,'K'), Tmax=(1074.38,'K')), NASAPolynomial(coeffs=[15.271,0.021725,-8.66327e-06,1.60999e-09,-1.13268e-13,63270.2,-48.0824], Tmin=(1074.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(Neopentyl) + radical(CCsJOH)"""),
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
    label = '[CH2][C]=CC(C)=[C]O(28972)',
    structure = SMILES('[CH2][C]=CC(C)=[C]O'),
    E0 = (407.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0109188,0.0727301,-7.08355e-05,3.47447e-08,-6.49901e-12,49125.8,29.7441], Tmin=(100,'K'), Tmax=(1465.77,'K')), NASAPolynomial(coeffs=[18.5133,0.014454,-3.23249e-06,3.74163e-10,-1.88615e-14,44537.9,-63.7416], Tmin=(1465.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([C]=[C]C)=CO(28973)',
    structure = SMILES('[CH2]C([C]=[C]C)=CO'),
    E0 = (399.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.17859,'amu*angstrom^2'), symmetry=1, barrier=(27.0982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17849,'amu*angstrom^2'), symmetry=1, barrier=(27.0959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18026,'amu*angstrom^2'), symmetry=1, barrier=(27.1364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17847,'amu*angstrom^2'), symmetry=1, barrier=(27.0954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.664516,0.0852975,-9.30673e-05,4.90825e-08,-9.65485e-12,48273.5,27.3517], Tmin=(100,'K'), Tmax=(1430.1,'K')), NASAPolynomial(coeffs=[22.0491,0.0093501,-3.83436e-07,-1.95292e-10,2.09395e-14,43046.7,-85.9143], Tmin=(1430.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CC([CH2])=[C]O(28974)',
    structure = SMILES('[CH2]C=CC([CH2])=[C]O'),
    E0 = (320.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,216.791],'cm^-1')),
        HinderedRotor(inertia=(0.85815,'amu*angstrom^2'), symmetry=1, barrier=(28.6202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858147,'amu*angstrom^2'), symmetry=1, barrier=(28.6202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858147,'amu*angstrom^2'), symmetry=1, barrier=(28.6202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858147,'amu*angstrom^2'), symmetry=1, barrier=(28.6202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646784,0.0585795,-1.30001e-05,-3.95423e-08,2.40044e-11,38719.6,27.2652], Tmin=(100,'K'), Tmax=(913.913,'K')), NASAPolynomial(coeffs=[20.6553,0.0112272,-1.29542e-06,7.4568e-11,-5.55411e-15,33382.7,-76.6529], Tmin=(913.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O)C=[C]C(28975)',
    structure = SMILES('[CH2]C(=[C]O)C=[C]C'),
    E0 = (440.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.862332,'amu*angstrom^2'), symmetry=1, barrier=(19.8267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862373,'amu*angstrom^2'), symmetry=1, barrier=(19.8277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862312,'amu*angstrom^2'), symmetry=1, barrier=(19.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861385,'amu*angstrom^2'), symmetry=1, barrier=(19.8049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.150027,0.0777747,-8.09453e-05,4.15852e-08,-8.12683e-12,53152.6,29.0413], Tmin=(100,'K'), Tmax=(1378.84,'K')), NASAPolynomial(coeffs=[20.0676,0.0125261,-2.78612e-06,3.25244e-10,-1.67233e-14,48204.3,-72.7178], Tmin=(1378.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=[C]C)=C[O](20119)',
    structure = SMILES('[CH2]C(C=[C]C)=C[O]'),
    E0 = (342.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07763,'amu*angstrom^2'), symmetry=1, barrier=(24.7769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07539,'amu*angstrom^2'), symmetry=1, barrier=(24.7254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07962,'amu*angstrom^2'), symmetry=1, barrier=(24.8225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530149,0.0634467,-3.12747e-05,-1.58434e-08,1.3841e-11,41307.8,24.8795], Tmin=(100,'K'), Tmax=(942.203,'K')), NASAPolynomial(coeffs=[19.5731,0.0144296,-3.90834e-06,6.45033e-10,-4.67862e-14,36306.6,-73.355], Tmin=(942.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]=C[C]1CC1O(28976)',
    structure = SMILES('[CH2][C]=C[C]1CC1O'),
    E0 = (447.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45922,0.040247,2.08753e-05,-5.90944e-08,2.62973e-11,53956.4,25.8636], Tmin=(100,'K'), Tmax=(962.886,'K')), NASAPolynomial(coeffs=[15.6745,0.019514,-6.52176e-06,1.20506e-09,-8.953e-14,49442.5,-51.4036], Tmin=(962.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(CCJ(C)CO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]1[CH]C(=C)C1O(28977)',
    structure = SMILES('[CH2][C]1[CH]C(=C)C1O'),
    E0 = (393.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.6806,-0.0203722,0.000139537,-1.31459e-07,3.03686e-11,46975.1,-14.7349], Tmin=(100,'K'), Tmax=(1695.76,'K')), NASAPolynomial(coeffs=[70.1873,0.0287694,-7.15609e-05,1.74324e-08,-1.29754e-12,-454.247,-413.588], Tmin=(1695.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Isobutyl) + radical(Allyl_S) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C=C1[CH][C]CC1O(28978)',
    structure = SMILES('C=C1[CH][C]CC1O'),
    E0 = (409.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80626,0.0263593,6.7276e-05,-1.08406e-07,4.35225e-11,49369.1,21.6161], Tmin=(100,'K'), Tmax=(965.395,'K')), NASAPolynomial(coeffs=[16.7916,0.0190417,-6.45768e-06,1.28135e-09,-1.01212e-13,43923.4,-63.3708], Tmin=(965.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclopentane) + radical(Allyl_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C[C]([CH2])C[O](20127)',
    structure = SMILES('[CH2][C]=C[C]([CH2])C[O]'),
    E0 = (728.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,360,370,350,1430.53,1430.61,1430.89],'cm^-1')),
        HinderedRotor(inertia=(0.0431103,'amu*angstrom^2'), symmetry=1, barrier=(62.6298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271403,'amu*angstrom^2'), symmetry=1, barrier=(6.24009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271495,'amu*angstrom^2'), symmetry=1, barrier=(6.24221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.72373,'amu*angstrom^2'), symmetry=1, barrier=(62.6238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45939,0.0594896,-5.75807e-05,3.6908e-08,-1.06298e-11,87715.6,30.0029], Tmin=(100,'K'), Tmax=(812.172,'K')), NASAPolynomial(coeffs=[5.99199,0.0371661,-1.63511e-05,3.0647e-09,-2.1223e-13,86979.4,9.07872], Tmin=(812.172,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(728.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P) + radical(CCJ(C)CO) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([CH2])[CH]O(20129)',
    structure = SMILES('[CH2][C]=[C]C([CH2])[CH]O'),
    E0 = (768.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,1670,1700,300,440,1380,1390,370,380,2900,435,314.671,317.23],'cm^-1')),
        HinderedRotor(inertia=(0.00146418,'amu*angstrom^2'), symmetry=1, barrier=(9.77,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00165686,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137672,'amu*angstrom^2'), symmetry=1, barrier=(9.76705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02139,'amu*angstrom^2'), symmetry=1, barrier=(73.1129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00644,'amu*angstrom^2'), symmetry=1, barrier=(73.1964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375529,0.0841306,-0.000120075,9.38943e-08,-2.88974e-11,92547.8,32.1559], Tmin=(100,'K'), Tmax=(889.029,'K')), NASAPolynomial(coeffs=[11.0391,0.0279329,-1.13886e-05,1.99355e-09,-1.3011e-13,90976.6,-16.2088], Tmin=(889.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(768.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Isobutyl) + radical(Cds_S) + radical(Allyl_P) + radical(CCsJOH)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3932.86,'J/mol'), sigma=(6.48948,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=614.30 K, Pc=32.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479077,0.0662053,-4.26455e-05,-3.6233e-09,9.60913e-12,56799,28.354], Tmin=(100,'K'), Tmax=(936.92,'K')), NASAPolynomial(coeffs=[19.3344,0.0137671,-3.61804e-06,5.73436e-10,-4.04467e-14,52034.2,-67.9562], Tmin=(936.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
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
    label = '[CH2]C(=[C]O)C[C]=C(28593)',
    structure = SMILES('[CH2]C(=[C]O)C[C]=C'),
    E0 = (469.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,372.331,372.335],'cm^-1')),
        HinderedRotor(inertia=(0.171824,'amu*angstrom^2'), symmetry=1, barrier=(16.9032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171824,'amu*angstrom^2'), symmetry=1, barrier=(16.9032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171821,'amu*angstrom^2'), symmetry=1, barrier=(16.9032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171824,'amu*angstrom^2'), symmetry=1, barrier=(16.9032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.539365,0.0676161,-5.40979e-05,1.36622e-08,2.1004e-12,56658.7,29.7431], Tmin=(100,'K'), Tmax=(958.427,'K')), NASAPolynomial(coeffs=[17.2122,0.0173346,-5.61404e-06,9.51028e-10,-6.52013e-14,52576.2,-54.6106], Tmin=(958.427,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC([CH2])=CO(19986)',
    structure = SMILES('[CH]=[C]CC([CH2])=CO'),
    E0 = (477.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.958427,'amu*angstrom^2'), symmetry=1, barrier=(22.0361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.956433,'amu*angstrom^2'), symmetry=1, barrier=(21.9903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.95574,'amu*angstrom^2'), symmetry=1, barrier=(21.9744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955606,'amu*angstrom^2'), symmetry=1, barrier=(21.9713,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3874.43,'J/mol'), sigma=(6.44487,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.18 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.694207,0.0803134,-8.09633e-05,3.89328e-08,-6.96019e-12,57599.4,31.4631], Tmin=(100,'K'), Tmax=(1590.54,'K')), NASAPolynomial(coeffs=[22.5117,0.00849225,-5.35427e-07,-9.9197e-11,1.1194e-14,51920.1,-85.9085], Tmin=(1590.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=C[O])C[C]=C(18133)',
    structure = SMILES('[CH2]C(=C[O])C[C]=C'),
    E0 = (371.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,353.484,353.531,353.538],'cm^-1')),
        HinderedRotor(inertia=(0.25096,'amu*angstrom^2'), symmetry=1, barrier=(22.2534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250923,'amu*angstrom^2'), symmetry=1, barrier=(22.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250941,'amu*angstrom^2'), symmetry=1, barrier=(22.2536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.748913,0.0587906,-2.34521e-05,-1.93988e-08,1.37989e-11,44834.4,27.2727], Tmin=(100,'K'), Tmax=(959.674,'K')), NASAPolynomial(coeffs=[18.0392,0.0169982,-5.44977e-06,9.68084e-10,-7.02385e-14,40121.7,-62.6948], Tmin=(959.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C1CC1=CO(28822)',
    structure = SMILES('C=[C]C1CC1=CO'),
    E0 = (287.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.645602,0.0597939,-2.12846e-05,-2.78895e-08,1.88585e-11,34735.4,22.1879], Tmin=(100,'K'), Tmax=(927.364,'K')), NASAPolynomial(coeffs=[20.045,0.0120966,-2.3291e-06,3.18615e-10,-2.37844e-14,29590.3,-78.2813], Tmin=(927.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=C[C]=C)CO(28979)',
    structure = SMILES('[CH]C(=C[C]=C)CO'),
    E0 = (435.385,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.668872,0.0722785,-6.86825e-05,3.72317e-08,-8.36231e-12,52485.5,26.9764], Tmin=(100,'K'), Tmax=(1064.65,'K')), NASAPolynomial(coeffs=[11.2179,0.0326449,-1.28423e-05,2.26553e-09,-1.51609e-13,50239.3,-24.5774], Tmin=(1064.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C(C=C[CH2])=CO(28980)',
    structure = SMILES('[CH]C(C=C[CH2])=CO'),
    E0 = (300.253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06407,'amu*angstrom^2'), symmetry=1, barrier=(47.4571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06563,'amu*angstrom^2'), symmetry=1, barrier=(47.4929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06415,'amu*angstrom^2'), symmetry=1, barrier=(47.4588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06223,'amu*angstrom^2'), symmetry=1, barrier=(47.4148,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440458,0.0591975,5.88605e-06,-6.58021e-08,3.46055e-11,36258,25.7056], Tmin=(100,'K'), Tmax=(912.715,'K')), NASAPolynomial(coeffs=[22.4076,0.0133102,-1.50449e-06,7.79e-11,-6.05088e-15,30149.4,-89.7634], Tmin=(912.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C(C=[C]C)=CO(28981)',
    structure = SMILES('[CH]C(C=[C]C)=CO'),
    E0 = (420.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.95078,'amu*angstrom^2'), symmetry=1, barrier=(44.8523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95315,'amu*angstrom^2'), symmetry=1, barrier=(44.9069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95884,'amu*angstrom^2'), symmetry=1, barrier=(45.0376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95738,'amu*angstrom^2'), symmetry=1, barrier=(45.004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.833872,0.0838102,-8.01244e-05,3.77002e-08,-6.69049e-12,50712.3,29.2106], Tmin=(100,'K'), Tmax=(1583.61,'K')), NASAPolynomial(coeffs=[21.8223,0.0143633,-2.76905e-06,2.62328e-10,-1.10027e-14,45068.9,-85.6694], Tmin=(1583.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]1C=C([CH]O)C1(28982)',
    structure = SMILES('[CH2][C]1C=C([CH]O)C1'),
    E0 = (374.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76694,0.0295785,5.60487e-05,-1.00049e-07,4.25789e-11,45186.4,25.1522], Tmin=(100,'K'), Tmax=(929.38,'K')), NASAPolynomial(coeffs=[16.4108,0.0165227,-3.5308e-06,5.42246e-10,-4.19805e-14,40306.4,-56.0338], Tmin=(929.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=CCJO) + radical(Allyl_T) + radical(Isobutyl)"""),
)

species(
    label = 'OC=C1[CH][C]CC1(28983)',
    structure = SMILES('OC=C1[CH][C]CC1'),
    E0 = (373.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4496,0.0277522,8.48994e-05,-1.43378e-07,6.05234e-11,45060.1,20.1812], Tmin=(100,'K'), Tmax=(934.345,'K')), NASAPolynomial(coeffs=[22.9593,0.00806209,2.87454e-07,-7.98451e-11,-6.70488e-15,37880.6,-99.0405], Tmin=(934.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(methylenecyclopentane) + radical(Allyl_S) + radical(CCJ2_triplet)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.5,'J/mol'), sigma=(6.3848,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.16 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00748,0.0534506,-1.61976e-05,-2.04947e-08,1.22809e-11,52644,30.7954], Tmin=(100,'K'), Tmax=(999.958,'K')), NASAPolynomial(coeffs=[16.1987,0.020262,-7.78247e-06,1.47611e-09,-1.07528e-13,48227.1,-49.3871], Tmin=(999.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=CO)C[C]=C(28598)',
    structure = SMILES('[CH]C(=CO)C[C]=C'),
    E0 = (449.428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.344897,0.0680946,-3.47321e-05,-1.31942e-08,1.29382e-11,54196.6,28.141], Tmin=(100,'K'), Tmax=(938.449,'K')), NASAPolynomial(coeffs=[18.8998,0.0195276,-5.88661e-06,9.69367e-10,-6.69432e-14,49370.1,-67.3571], Tmin=(938.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(=C)C=O(18125)',
    structure = SMILES('C=[C]CC(=C)C=O'),
    E0 = (195.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2782.5,750,1395,475,1775,1000,397.32,397.339],'cm^-1')),
        HinderedRotor(inertia=(0.11068,'amu*angstrom^2'), symmetry=1, barrier=(12.3977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110703,'amu*angstrom^2'), symmetry=1, barrier=(12.3974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110694,'amu*angstrom^2'), symmetry=1, barrier=(12.3969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78024,0.0535918,-3.60573e-05,1.17702e-08,-1.60459e-12,23622.5,23.8722], Tmin=(100,'K'), Tmax=(1580.6,'K')), NASAPolynomial(coeffs=[10.5262,0.0314586,-1.50527e-05,2.91091e-09,-2.03338e-13,20857.7,-22.3258], Tmin=(1580.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1C(=C)C1O(28827)',
    structure = SMILES('C=[C]C1C(=C)C1O'),
    E0 = (323.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.929445,0.0592335,-4.1662e-05,1.03564e-08,6.2036e-13,39047.5,23.1922], Tmin=(100,'K'), Tmax=(1065.75,'K')), NASAPolynomial(coeffs=[14.3928,0.022218,-8.58685e-06,1.56597e-09,-1.08928e-13,35410.2,-46.2194], Tmin=(1065.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C#CC([CH2])[CH]O(28984)',
    structure = SMILES('[CH2]C#CC([CH2])[CH]O'),
    E0 = (456.455,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,2100,2250,500,550,1380,1390,370,380,2900,435,373.397,373.397],'cm^-1')),
        HinderedRotor(inertia=(0.0180902,'amu*angstrom^2'), symmetry=1, barrier=(71.6645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11281,'amu*angstrom^2'), symmetry=1, barrier=(11.1613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724328,'amu*angstrom^2'), symmetry=1, barrier=(71.6645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112811,'amu*angstrom^2'), symmetry=1, barrier=(11.1613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180902,'amu*angstrom^2'), symmetry=1, barrier=(71.6645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.92728,0.0701908,-7.52662e-05,3.80869e-08,-4.1295e-12,55007.3,29.4516], Tmin=(100,'K'), Tmax=(714.722,'K')), NASAPolynomial(coeffs=[11.1241,0.0259176,-9.20017e-06,1.50875e-09,-9.56278e-14,53223,-18.6035], Tmin=(714.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = 'C#C[CH]C([CH2])[CH]O(28985)',
    structure = SMILES('C#C[CH]C([CH2])[CH]O'),
    E0 = (464.743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3000,3050,390,425,1340,1360,335,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2175,525,3615,1277.5,1000,244.519,2884.48],'cm^-1')),
        HinderedRotor(inertia=(0.288159,'amu*angstrom^2'), symmetry=1, barrier=(12.2601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75169,'amu*angstrom^2'), symmetry=1, barrier=(74.5696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73833,'amu*angstrom^2'), symmetry=1, barrier=(74.5936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126423,'amu*angstrom^2'), symmetry=1, barrier=(74.5791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74866,'amu*angstrom^2'), symmetry=1, barrier=(74.6247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.289275,0.082837,-0.000112212,8.1027e-08,-2.24359e-11,56028.1,29.3172], Tmin=(100,'K'), Tmax=(1017.94,'K')), NASAPolynomial(coeffs=[12.6584,0.0227712,-6.81202e-06,9.38216e-10,-5.0116e-14,54103.7,-27.66], Tmin=(1017.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOH) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2][CH][C]=C([CH2])[CH]O(28986)',
    structure = SMILES('[CH2][CH][C]=C([CH2])[CH]O'),
    E0 = (608.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,1685,370,350,440,435,1725,425.961,430.455],'cm^-1')),
        HinderedRotor(inertia=(0.00091041,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000909611,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000907774,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00090502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.377331,'amu*angstrom^2'), symmetry=1, barrier=(48.0622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908754,0.0582715,-3.3175e-05,-6.00734e-10,4.82114e-12,73276.2,30.2165], Tmin=(100,'K'), Tmax=(1027.63,'K')), NASAPolynomial(coeffs=[15.0178,0.0224378,-8.7272e-06,1.61144e-09,-1.13688e-13,69368.7,-43.1389], Tmin=(1027.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(608.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(RCCJ) + radical(Allyl_S) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH][CH]C=C([CH2])[CH]O(28987)',
    structure = SMILES('[CH][CH]C=C([CH2])[CH]O'),
    E0 = (613.384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.875955,0.0557298,-1.66227e-05,-2.33571e-08,1.41532e-11,73897.1,28.384], Tmin=(100,'K'), Tmax=(979.711,'K')), NASAPolynomial(coeffs=[17.0609,0.0195679,-7.06429e-06,1.3095e-09,-9.51973e-14,69290,-56.6954], Tmin=(979.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(Allyl_P) + radical(CCJ2_triplet) + radical(Allyl_S)"""),
)

species(
    label = '[CH2][C][CH]O(28660)',
    structure = SMILES('[CH2][C][CH]O'),
    E0 = (553.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,180,2262.1],'cm^-1')),
        HinderedRotor(inertia=(0.0475559,'amu*angstrom^2'), symmetry=1, barrier=(11.7885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.515056,'amu*angstrom^2'), symmetry=1, barrier=(11.8422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400962,'amu*angstrom^2'), symmetry=1, barrier=(46.0504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12807,0.0437773,-6.51486e-05,5.20447e-08,-1.62809e-11,66624.5,18.0235], Tmin=(100,'K'), Tmax=(877.804,'K')), NASAPolynomial(coeffs=[7.75761,0.0135713,-5.75184e-06,1.02548e-09,-6.75894e-14,65811.6,-7.40289], Tmin=(877.804,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCsJOH) + radical(RCCJ)"""),
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
    E0 = (302.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (513.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (515.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (502.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (617.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (654.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (520.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (514.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (501.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (529.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (522.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (405.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (451.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (493.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (445.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (672.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (710.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (758.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (504.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (492.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (489.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (363.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (378.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (325.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (339.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (307.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (307.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (308.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (740.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (853.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (819.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (915.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (770.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (749.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (558.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (620.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (537.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (555.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (592.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (724.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (473.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (435.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (801.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (447.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (438.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (416.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (756.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (791.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (641.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (443.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (620.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (521.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (485.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (305.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (580.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (703.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (453.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (438.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (416.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (531.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (641.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (325.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (323.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (634.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (779.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (631.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (636.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (916.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=C=CO(12571)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(2.36501,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 2.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH2]C(C=C=C)=C[O](20110)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2834 used for Od_CO-CdH;HJ
Exact match found for rate rule [Od_CO-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C#CC([CH2])=CO(28959)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC=C([CH2])[CH]O(27790)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.156e+09,'cm^3/(mol*s)'), n=1.502, Ea=(9.92026,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 194 used for Ct-H_Ct-Cd;HJ
Exact match found for rate rule [Ct-H_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(28)', 'C=[C]C=C=CO(28960)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(10.4286,'m^3/(mol*s)'), n=1.99344, Ea=(11.9048,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds;YJ] for rate rule [Ca_Cds;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]O(5471)', 'C=[C]C=C=C(19283)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.40891,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=C[C]=C)C[O](20752)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2]C=[C]C([CH2])=CO(28961)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]C=CC([CH2])=CO(28744)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=[C][C]=C)CO(28962)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=[C]C(C)=CO(28963)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=CC(C)=C[O](20112)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.3e+09,'s^-1','*|/',2.51189), n=0.75, Ea=(97.0688,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;Cs_H_out_2H] for rate rule [R4HJ_1;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2]C=CC([CH2])=C[O](20116)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]C=C([CH2])CO(28756)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.1728e+06,'s^-1'), n=1.70245, Ea=(63.8935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cs_H_out_1H] + [R5Hall;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#C[CH][C](C)[CH]O(28745)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH2][C]=CC([CH2])=C[O](20121)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH2]C(=[C][C]=C)[CH]O(28964)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2][C]1C(O)C1[C]=C(28965)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.00692e+11,'s^-1'), n=0.347401, Ea=(204.477,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 202.4 to 204.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=[C]C1C[C]1[CH]O(28966)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.03537e+13,'s^-1'), n=-0.296394, Ea=(192.336,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 191.4 to 192.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2]C1([CH]O)[CH]C1=C(28801)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(189.653,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2]C#CC(=C)CO(28967)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2]C#CC(C)=CO(28968)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=[C]C=C(C)C=O(20136)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2]C=CC(=C)C=O(20126)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=C=C[C]1CC1O(28969)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=C1[CH]C(=C)C1O(28783)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_H/NonDeO;Ypri_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=C1[CH]C(=CO)C1(28772)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['OH(D)(132)', '[CH]C([CH2])=C[C]=C(18814)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad;Birad] for rate rule [O_pri_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]O(5471)', '[CH2][C]=C[C]=C(19278)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(T)(28)', '[CH2][C]=C[C]=CO(28970)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]=C(584)', '[CH]C([CH2])=CO(19098)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(8)', '[CH2][C]=CC([CH2])=[C]O(28971)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(8)', '[CH]C(C=[C][CH2])=CO(24296)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2]C1([CH]O)C=[C]C1(28743)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.17094e+10,'s^-1'), n=0.2405, Ea=(258.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5;doublebond_intra;radadd_intra_cs2H] + [R5;doublebond_intra_HNd;radadd_intra_cs] for rate rule [R5;doublebond_intra_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C=CO(12571)', '[CH][C]=C(18825)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][C]=CO(18753)', 'C#C[CH2](17441)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2][C]=CC(C)=[C]O(28972)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([C]=[C]C)=CO(28973)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C=CC([CH2])=[C]O(28974)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(=[C]O)C=[C]C(28975)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(C=[C]C)=C[O](20119)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_2H] for rate rule [R6HJ_4;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][C]=CO(18753)', '[CH][C]=C(18825)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2][C]=C[C]1CC1O(28976)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.03537e+13,'s^-1'), n=-0.296394, Ea=(147.69,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 145.5 to 147.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2][C]1[CH]C(=C)C1O(28977)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=C1[CH][C]CC1O(28978)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.4946e+10,'s^-1'), n=0.314866, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2][C]=C[C]([CH2])C[O](20127)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][C]=[C]C([CH2])[CH]O(20129)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=CC[C]=CO(20241)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2][C]=CC([CH2])C=O(18134)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C(=[C]O)C[C]=C(28593)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(=C[O])C[C]=C(18133)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=[C]C1CC1=CO(28822)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]C(=C[C]=C)CO(28979)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]C(C=C[CH2])=CO(28980)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH]C(C=[C]C)=CO(28981)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['[CH2][C]1C=C([CH]O)C1(28982)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['OC=C1[CH][C]CC1(28983)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(4.4946e+10,'s^-1'), n=0.314866, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH]C(=CO)C[C]=C(28598)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=[C]CC(=C)C=O(18125)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    products = ['C=[C]C1C(=C)C1O(28827)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(23.7316,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C#CC([CH2])[CH]O(28984)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction65',
    reactants = ['C#C[CH]C([CH2])[CH]O(28985)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_MMS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2][CH][C]=C([CH2])[CH]O(28986)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH][CH]C=C([CH2])[CH]O(28987)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction68',
    reactants = ['C#C[CH2](17441)', '[CH2][C][CH]O(28660)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '5229',
    isomers = [
        '[CH2]C([CH]O)=C[C]=C(20132)',
    ],
    reactants = [
        ('C=C=CO(12571)', 'C#C[CH2](17441)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5229',
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

