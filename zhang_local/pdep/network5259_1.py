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
    label = '[O]C=[C]C1CC1=CO(29139)',
    structure = SMILES('[O]C=[C]C1CC1=CO'),
    E0 = (220.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0215818,0.0685472,-2.21542e-05,-4.33127e-08,2.83969e-11,26664.8,25.6982], Tmin=(100,'K'), Tmax=(912.308,'K')), NASAPolynomial(coeffs=[26.5571,0.002819,2.69274e-06,-6.55026e-10,4.24508e-14,19716.7,-111.429], Tmin=(912.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Methylene_cyclopropane) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1([CH]O)C=C=CO1(29192)',
    structure = SMILES('[CH2]C1([CH]O)C=C=CO1'),
    E0 = (436.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31047,0.0936764,-0.000104103,5.37179e-08,-1.02091e-11,52659.8,28.7361], Tmin=(100,'K'), Tmax=(1497.1,'K')), NASAPolynomial(coeffs=[26.2762,0.00399859,1.75056e-06,-5.44644e-10,4.20515e-14,46189.6,-109.509], Tmin=(1497.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(CJC(C)OC) + radical(CCsJOH)"""),
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
    label = '[CH2]C(C=C=C=O)=CO(29193)',
    structure = SMILES('[CH2]C(C=C=C=O)=CO'),
    E0 = (116.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,350,440,435,1725,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.113415,'amu*angstrom^2'), symmetry=1, barrier=(2.60764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113358,'amu*angstrom^2'), symmetry=1, barrier=(2.60631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113482,'amu*angstrom^2'), symmetry=1, barrier=(2.60917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849024,0.0506011,-1.19058e-06,-5.69154e-08,3.20616e-11,14128.5,9.37802], Tmin=(100,'K'), Tmax=(909.314,'K')), NASAPolynomial(coeffs=[24.503,-0.00338641,5.28092e-06,-1.11186e-09,7.27562e-14,7756.97,-113.871], Tmin=(909.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Allyl_P)"""),
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
    label = '[CH2]C([CH]O)=CC#CO(29194)',
    structure = SMILES('[CH2]C([CH]O)=CC#CO'),
    E0 = (139.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,350,440,435,1725,2100,2250,500,550,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,832.1],'cm^-1')),
        HinderedRotor(inertia=(0.978564,'amu*angstrom^2'), symmetry=1, barrier=(22.4991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0458316,'amu*angstrom^2'), symmetry=1, barrier=(22.4989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243669,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0458244,'amu*angstrom^2'), symmetry=1, barrier=(22.5002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60734,0.0647603,-4.50787e-05,7.73936e-09,2.74013e-12,16937.2,28.9696], Tmin=(100,'K'), Tmax=(1023.31,'K')), NASAPolynomial(coeffs=[16.6858,0.0206709,-7.94905e-06,1.46455e-09,-1.03507e-13,12664.4,-53.7692], Tmin=(1023.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CtH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(CTCC=CCJ) + radical(C=CCJO)"""),
)

species(
    label = 'CC(=C=C=C[O])[CH]O(29195)',
    structure = SMILES('CC(=C=C=C[O])[CH]O'),
    E0 = (145.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,540,563.333,586.667,610,1970,2140,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18589,'amu*angstrom^2'), symmetry=1, barrier=(27.266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18558,'amu*angstrom^2'), symmetry=1, barrier=(27.2588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18541,'amu*angstrom^2'), symmetry=1, barrier=(27.255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.183499,0.0722122,-5.21433e-05,5.12969e-09,6.1149e-12,17594.1,28.9612], Tmin=(100,'K'), Tmax=(968.908,'K')), NASAPolynomial(coeffs=[20.369,0.0153339,-5.04373e-06,9.02226e-10,-6.53545e-14,12440.7,-74.1928], Tmin=(968.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(C=COJ)"""),
)

species(
    label = 'CC(=[C]O)C=C=C[O](29196)',
    structure = SMILES('CC(=[C]O)C=C=C[O]'),
    E0 = (160.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30933,0.0911281,-9.89568e-05,4.99445e-08,-9.2523e-12,19520.6,32.7798], Tmin=(100,'K'), Tmax=(1553.42,'K')), NASAPolynomial(coeffs=[25.6747,0.00416751,1.8901e-06,-5.77888e-10,4.41751e-14,13245.9,-102.5], Tmin=(1553.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=COJ)"""),
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
    label = '[CH2]C(=C=C=CO)[CH]O(29197)',
    structure = SMILES('[CH2]C(=C=C=CO)[CH]O'),
    E0 = (155.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,540,563.333,586.667,610,1970,2140,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.31275,'amu*angstrom^2'), symmetry=1, barrier=(30.1827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31351,'amu*angstrom^2'), symmetry=1, barrier=(30.2001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31482,'amu*angstrom^2'), symmetry=1, barrier=(30.2302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31299,'amu*angstrom^2'), symmetry=1, barrier=(30.1883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47095,0.0909629,-9.53036e-05,4.61731e-08,-8.20673e-12,18875.5,33.7951], Tmin=(100,'K'), Tmax=(1625.39,'K')), NASAPolynomial(coeffs=[26.2465,0.00441216,1.49498e-06,-4.71407e-10,3.53531e-14,12287.7,-105.937], Tmin=(1625.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = 'CC([CH][C]=C=O)=CO(29198)',
    structure = SMILES('CC([CH][C]=C=O)=CO'),
    E0 = (104.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2120,512.5,787.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.998249,'amu*angstrom^2'), symmetry=1, barrier=(22.9517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999111,'amu*angstrom^2'), symmetry=1, barrier=(22.9715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99998,'amu*angstrom^2'), symmetry=1, barrier=(22.9915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999486,'amu*angstrom^2'), symmetry=1, barrier=(22.9802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.260037,0.0891512,-0.000105891,6.25687e-08,-1.43007e-11,12665.8,28.4371], Tmin=(100,'K'), Tmax=(1080.79,'K')), NASAPolynomial(coeffs=[19.1585,0.017283,-6.14667e-06,1.04269e-09,-6.89283e-14,8468.37,-66.7545], Tmin=(1080.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCJC(C)=C=O) + radical(CCCJ=C=O)"""),
)

species(
    label = '[CH2]C(=[C]O)C=C=CO(29199)',
    structure = SMILES('[CH2]C(=[C]O)C=C=CO'),
    E0 = (170.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,540,610,2055,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.20418,'amu*angstrom^2'), symmetry=1, barrier=(27.6864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20295,'amu*angstrom^2'), symmetry=1, barrier=(27.6582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2019,'amu*angstrom^2'), symmetry=1, barrier=(27.6341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20536,'amu*angstrom^2'), symmetry=1, barrier=(27.7136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.28371,0.10217,-0.000116541,5.96789e-08,-1.09776e-11,20771.8,35.1523], Tmin=(100,'K'), Tmax=(1614.34,'K')), NASAPolynomial(coeffs=[29.0051,-0.00258366,6.09136e-06,-1.41121e-09,1.00812e-13,14217.4,-119.795], Tmin=(1614.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(C=C=CO)=C[O](25379)',
    structure = SMILES('[CH2]C(C=C=CO)=C[O]'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.42989,0.0741445,-1.93668e-05,-5.92742e-08,3.73878e-11,8875.35,26.7673], Tmin=(100,'K'), Tmax=(904.884,'K')), NASAPolynomial(coeffs=[31.7552,-0.00493166,6.95647e-06,-1.48733e-09,9.9297e-14,463.293,-139.586], Tmin=(904.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.2858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Allyl_P)"""),
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
    label = '[CH]C(=C)C=C=C[O](22618)',
    structure = SMILES('[CH]C(=C)C=C=C[O]'),
    E0 = (500.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.00191,'amu*angstrom^2'), symmetry=1, barrier=(46.0278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00636,'amu*angstrom^2'), symmetry=1, barrier=(46.1302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3964.84,'J/mol'), sigma=(6.40498,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=619.30 K, Pc=34.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61675,0.0607309,-2.30026e-05,-2.49436e-08,1.72533e-11,60302.2,23.9773], Tmin=(100,'K'), Tmax=(936.072,'K')), NASAPolynomial(coeffs=[19.4917,0.0144273,-3.85094e-06,6.21069e-10,-4.4817e-14,55263.5,-73.8753], Tmin=(936.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(500.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
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
    label = '[CH2]C(=C=C=C[O])[CH]O(29200)',
    structure = SMILES('[CH2]C(=C=C=C[O])[CH]O'),
    E0 = (296.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,540,563.333,586.667,610,1970,2140,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.48324,'amu*angstrom^2'), symmetry=1, barrier=(34.1027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48144,'amu*angstrom^2'), symmetry=1, barrier=(34.0613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48516,'amu*angstrom^2'), symmetry=1, barrier=(34.1467,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202067,0.0694859,-4.27151e-05,-9.48979e-09,1.27297e-11,35816.9,29.1755], Tmin=(100,'K'), Tmax=(945.955,'K')), NASAPolynomial(coeffs=[22.5106,0.00953678,-2.17548e-06,3.64127e-10,-2.94325e-14,30057.9,-85.3423], Tmin=(945.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=CCJO) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(=[C]O)C=C=C[O](29201)',
    structure = SMILES('[CH2]C(=[C]O)C=C=C[O]'),
    E0 = (312.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30717,'amu*angstrom^2'), symmetry=1, barrier=(30.0545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31029,'amu*angstrom^2'), symmetry=1, barrier=(30.1261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30936,'amu*angstrom^2'), symmetry=1, barrier=(30.1048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.61803,0.0921588,-0.000102267,5.14043e-08,-9.35014e-12,37757.9,34.1761], Tmin=(100,'K'), Tmax=(1616.33,'K')), NASAPolynomial(coeffs=[26.7076,-0.000196568,4.09626e-06,-9.85413e-10,7.07641e-14,31508.5,-107.086], Tmin=(1616.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([CH][C]=C=O)=CO(29202)',
    structure = SMILES('[CH2]C([CH][C]=C=O)=CO'),
    E0 = (255.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24441,'amu*angstrom^2'), symmetry=1, barrier=(28.6115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24654,'amu*angstrom^2'), symmetry=1, barrier=(28.6605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24291,'amu*angstrom^2'), symmetry=1, barrier=(28.577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25135,'amu*angstrom^2'), symmetry=1, barrier=(28.7711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.24827,0.086508,-9.67793e-05,4.84238e-08,-7.928e-12,30888.9,28.6758], Tmin=(100,'K'), Tmax=(943.475,'K')), NASAPolynomial(coeffs=[21.2165,0.0116267,-3.35921e-06,5.23611e-10,-3.45812e-14,26121.1,-77.4328], Tmin=(943.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CCCJ=C=O) + radical(CCJC(C)=C=O)"""),
)

species(
    label = '[O]C=C=C[C]1CC1O(29203)',
    structure = SMILES('[O]C=C=C[C]1CC1O'),
    E0 = (167.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909524,0.046537,2.70873e-05,-8.16534e-08,3.8211e-11,20297.6,27.4056], Tmin=(100,'K'), Tmax=(938.38,'K')), NASAPolynomial(coeffs=[22.1929,0.0100876,-1.40589e-06,2.25867e-10,-2.39398e-14,13913.6,-86.6531], Tmin=(938.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(C=COJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C(=CO)C=C1[CH]O1(29204)',
    structure = SMILES('[CH2]C(=CO)C=C1[CH]O1'),
    E0 = (98.4679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0546438,0.0534841,5.60628e-05,-1.49533e-07,7.24733e-11,12022.4,24.002], Tmin=(100,'K'), Tmax=(906.795,'K')), NASAPolynomial(coeffs=[37.0827,-0.0154704,1.32049e-05,-2.6582e-09,1.74248e-13,1387.04,-173.036], Tmin=(906.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.4679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(methyleneoxirane) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = '[O]C=C1[CH]C(=CO)C1(29205)',
    structure = SMILES('[O]C=C1[CH]C(=CO)C1'),
    E0 = (21.7384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976784,0.0266811,0.000122179,-2.08619e-07,9.10607e-11,2760.59,25.3684], Tmin=(100,'K'), Tmax=(911.697,'K')), NASAPolynomial(coeffs=[33.0413,-0.00993347,1.12027e-05,-2.26824e-09,1.44566e-13,-7410.97,-150.079], Tmin=(911.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.7384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=COJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2]C(=CO)C1[C]=CO1(29206)',
    structure = SMILES('[CH2]C(=CO)C1[C]=CO1'),
    E0 = (204.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00208546,0.0566956,3.62222e-05,-1.1959e-07,5.90117e-11,24737.7,26.054], Tmin=(100,'K'), Tmax=(915.491,'K')), NASAPolynomial(coeffs=[33.7135,-0.00858036,8.79196e-06,-1.75668e-09,1.11204e-13,15128.1,-152.379], Tmin=(915.491,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[O]C1[C]=CC(=CO)C1(29207)',
    structure = SMILES('[O]C1[C]=CC(=CO)C1'),
    E0 = (203.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.718649,0.0481657,3.13597e-05,-9.31636e-08,4.41799e-11,24582.6,25.7281], Tmin=(100,'K'), Tmax=(929.344,'K')), NASAPolynomial(coeffs=[25.0342,0.00512374,1.38236e-06,-3.19447e-10,1.33335e-14,17402.3,-104.114], Tmin=(929.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(3-Methylenecyclopentene) + radical(cyclopentene-vinyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1=C[C]=COC1O(29208)',
    structure = SMILES('[CH2]C1=C[C]=COC1O'),
    E0 = (19.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759283,0.0505804,2.06597e-05,-7.77011e-08,3.7795e-11,2516.98,21.6667], Tmin=(100,'K'), Tmax=(925.704,'K')), NASAPolynomial(coeffs=[22.3601,0.0105652,-9.03736e-07,5.37654e-11,-8.69975e-15,-3766.91,-93.2177], Tmin=(925.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CC(C=C=C=O)=CO(29209)',
    structure = SMILES('CC(C=C=C=O)=CO'),
    E0 = (-35.1168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.845992,0.053152,-1.00577e-05,-4.29163e-08,2.56485e-11,-4094.89,9.10761], Tmin=(100,'K'), Tmax=(914.957,'K')), NASAPolynomial(coeffs=[22.2449,0.00260739,2.2998e-06,-5.47189e-10,3.46365e-14,-9810.87,-102.065], Tmin=(914.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.1168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C(=[C]C=C[O])[CH]O(29210)',
    structure = SMILES('[CH2]C(=[C]C=C[O])[CH]O'),
    E0 = (232.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,337.008,337.01,337.011],'cm^-1')),
        HinderedRotor(inertia=(0.516114,'amu*angstrom^2'), symmetry=1, barrier=(41.5969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516115,'amu*angstrom^2'), symmetry=1, barrier=(41.5969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516112,'amu*angstrom^2'), symmetry=1, barrier=(41.5968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516111,'amu*angstrom^2'), symmetry=1, barrier=(41.5969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.349528,0.0633278,-1.43925e-05,-4.40838e-08,2.69199e-11,28137.2,30.4381], Tmin=(100,'K'), Tmax=(912.266,'K')), NASAPolynomial(coeffs=[22.8858,0.00964977,-3.47719e-07,-1.12015e-10,6.94415e-15,22147.2,-86.511], Tmin=(912.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=CCJO) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C([CH]C=[C][O])=CO(29211)',
    structure = SMILES('[CH2]C([CH]C=[C][O])=CO'),
    E0 = (266.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,347.315,347.457,347.517],'cm^-1')),
        HinderedRotor(inertia=(0.349323,'amu*angstrom^2'), symmetry=1, barrier=(29.9013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348905,'amu*angstrom^2'), symmetry=1, barrier=(29.9048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348701,'amu*angstrom^2'), symmetry=1, barrier=(29.9069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348789,'amu*angstrom^2'), symmetry=1, barrier=(29.9069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281931,0.0591907,8.19021e-06,-7.68388e-08,4.09641e-11,32212.4,31.1355], Tmin=(100,'K'), Tmax=(909.091,'K')), NASAPolynomial(coeffs=[27.0596,0.00155593,3.97909e-06,-9.24372e-10,6.04345e-14,24856.7,-109.178], Tmin=(909.091,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CCJC=C) + radical(C=COJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O)C[C]=C[O](29212)',
    structure = SMILES('[CH2]C(=[C]O)C[C]=C[O]'),
    E0 = (402.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,224.721,224.721,224.725],'cm^-1')),
        HinderedRotor(inertia=(0.60082,'amu*angstrom^2'), symmetry=1, barrier=(21.5311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600836,'amu*angstrom^2'), symmetry=1, barrier=(21.5311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600847,'amu*angstrom^2'), symmetry=1, barrier=(21.5311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600836,'amu*angstrom^2'), symmetry=1, barrier=(21.5312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02339,0.0872686,-9.23939e-05,4.5991e-08,-8.47244e-12,48629.3,36.6338], Tmin=(100,'K'), Tmax=(1542.55,'K')), NASAPolynomial(coeffs=[24.5449,0.00635422,4.98376e-07,-2.96892e-10,2.47293e-14,42479.7,-92.1655], Tmin=(1542.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=[C]O)[CH]C=C[O](29213)',
    structure = SMILES('[CH2]C(=[C]O)[CH]C=C[O]'),
    E0 = (266.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,347.486,347.486,347.486],'cm^-1')),
        HinderedRotor(inertia=(0.349016,'amu*angstrom^2'), symmetry=1, barrier=(29.905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349013,'amu*angstrom^2'), symmetry=1, barrier=(29.905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349015,'amu*angstrom^2'), symmetry=1, barrier=(29.905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349015,'amu*angstrom^2'), symmetry=1, barrier=(29.905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281974,0.0591902,8.19226e-06,-7.68417e-08,4.09655e-11,32212.4,31.1354], Tmin=(100,'K'), Tmax=(909.087,'K')), NASAPolynomial(coeffs=[27.0595,0.00155614,3.97897e-06,-9.24342e-10,6.0432e-14,24856.7,-109.177], Tmin=(909.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CJO) + radical(C=COJ) + radical(C=CCJC=C)"""),
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
    label = '[CH2]C(=[C]O)C=[C]C[O](29214)',
    structure = SMILES('[CH2]C(=[C]O)C=[C]C[O]'),
    E0 = (513.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,268.569,268.576,268.58],'cm^-1')),
        HinderedRotor(inertia=(0.00233706,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289595,'amu*angstrom^2'), symmetry=1, barrier=(14.8229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289575,'amu*angstrom^2'), symmetry=1, barrier=(14.8229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.451963,'amu*angstrom^2'), symmetry=1, barrier=(23.1346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.256959,0.0825455,-9.91806e-05,6.23138e-08,-1.54773e-11,61854.4,32.1991], Tmin=(100,'K'), Tmax=(986.226,'K')), NASAPolynomial(coeffs=[14.7879,0.0236091,-9.54018e-06,1.71808e-09,-1.16572e-13,58988.2,-37.7026], Tmin=(986.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(Allyl_P) + radical(Cds_S) + radical(CCOJ)"""),
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
    label = '[CH2]C(=C[C]=[C][O])CO(29215)',
    structure = SMILES('[CH2]C(=C[C]=[C][O])CO'),
    E0 = (355.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,252,252.223,2092.47],'cm^-1')),
        HinderedRotor(inertia=(1.93744,'amu*angstrom^2'), symmetry=1, barrier=(87.0398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93136,'amu*angstrom^2'), symmetry=1, barrier=(87.0144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92885,'amu*angstrom^2'), symmetry=1, barrier=(86.9515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372825,'amu*angstrom^2'), symmetry=1, barrier=(16.8771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330068,0.0764482,-8.26044e-05,4.7472e-08,-1.07672e-11,42852.9,33.0358], Tmin=(100,'K'), Tmax=(1081.64,'K')), NASAPolynomial(coeffs=[14.9589,0.0223472,-7.57512e-06,1.22594e-09,-7.79055e-14,39688.4,-38.6875], Tmin=(1081.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(C=CJO) + radical(C=COJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4556.79,'J/mol'), sigma=(7.12408,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=711.76 K, Pc=28.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0688376,0.0724563,-3.61978e-05,-2.66638e-08,2.18221e-11,23140.1,29.8906], Tmin=(100,'K'), Tmax=(916.915,'K')), NASAPolynomial(coeffs=[25.9182,0.00423282,1.55867e-06,-4.19886e-10,2.63005e-14,16476.8,-103.576], Tmin=(916.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S)"""),
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
    label = '[O]C=C=C[C]=CO(29216)',
    structure = SMILES('[O]C=C=C[C]=CO'),
    E0 = (157.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.47341,'amu*angstrom^2'), symmetry=1, barrier=(33.8766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4862,'amu*angstrom^2'), symmetry=1, barrier=(34.1707,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37295,0.0830796,-9.52177e-05,4.84908e-08,-8.77246e-12,19151.4,28.8859], Tmin=(100,'K'), Tmax=(1673.22,'K')), NASAPolynomial(coeffs=[23.8011,-0.00362147,6.28216e-06,-1.42271e-09,1.00575e-13,14439.4,-94.4294], Tmin=(1673.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=COJ)"""),
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
    label = '[CH]C(C=C=C[O])=CO(29217)',
    structure = SMILES('[CH]C(C=C=C[O])=CO'),
    E0 = (291.471,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.89686,'amu*angstrom^2'), symmetry=1, barrier=(43.6126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89699,'amu*angstrom^2'), symmetry=1, barrier=(43.6156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89539,'amu*angstrom^2'), symmetry=1, barrier=(43.5788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.419004,0.0764389,-2.71151e-05,-4.69911e-08,3.18217e-11,35234.5,27.5565], Tmin=(100,'K'), Tmax=(904.926,'K')), NASAPolynomial(coeffs=[29.3358,0.00133077,3.8694e-06,-9.24903e-10,6.24627e-14,27539.4,-125.784], Tmin=(904.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C1C(O)C1[C]=C[O](29165)',
    structure = SMILES('C=C1C(O)C1[C]=C[O]'),
    E0 = (256.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.359748,0.0674071,-4.0847e-05,-6.63902e-09,1.05067e-11,30974.5,26.5042], Tmin=(100,'K'), Tmax=(955.419,'K')), NASAPolynomial(coeffs=[20.3419,0.0138727,-4.09241e-06,7.1516e-10,-5.27723e-14,25781.3,-76.182], Tmin=(955.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=COJ) + radical(Cds_S)"""),
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
    label = 'C=C(C=C=C[O])C[O](26530)',
    structure = SMILES('C=C(C=C=C[O])C[O]'),
    E0 = (202.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03833,'amu*angstrom^2'), symmetry=1, barrier=(23.8732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03722,'amu*angstrom^2'), symmetry=1, barrier=(23.8477,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283801,0.0742833,-7.19283e-05,3.53731e-08,-6.84906e-12,24452.9,28.6816], Tmin=(100,'K'), Tmax=(1258,'K')), NASAPolynomial(coeffs=[17.2036,0.0204844,-7.78003e-06,1.37825e-09,-9.33237e-14,20195.9,-56.83], Tmin=(1258,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=C=C=C[O])CO(29218)',
    structure = SMILES('[CH2]C(=C=C=C[O])CO'),
    E0 = (179.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,563.333,586.667,610,1970,2140,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10937,'amu*angstrom^2'), symmetry=1, barrier=(25.5066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11345,'amu*angstrom^2'), symmetry=1, barrier=(25.6005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11252,'amu*angstrom^2'), symmetry=1, barrier=(25.5791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0950561,0.0815708,-8.6154e-05,4.53097e-08,-9.24319e-12,21714.6,29.6456], Tmin=(100,'K'), Tmax=(1206.88,'K')), NASAPolynomial(coeffs=[19.3713,0.0170522,-5.96494e-06,1.01384e-09,-6.74243e-14,17015.9,-67.9288], Tmin=(1206.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C(C=C=C[O])CO(29219)',
    structure = SMILES('[CH]=C(C=C=C[O])CO'),
    E0 = (223.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,540,610,2055,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.99459,'amu*angstrom^2'), symmetry=1, barrier=(22.8676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989034,'amu*angstrom^2'), symmetry=1, barrier=(22.7398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990489,'amu*angstrom^2'), symmetry=1, barrier=(22.7733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.510009,0.0847857,-9.12435e-05,4.74481e-08,-9.34388e-12,27060.3,31.0728], Tmin=(100,'K'), Tmax=(1360.94,'K')), NASAPolynomial(coeffs=[22.4297,0.0106977,-2.23932e-06,2.50248e-10,-1.27619e-14,21433.6,-84.3999], Tmin=(1360.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(=C[C]=C=O)CO(29220)',
    structure = SMILES('[CH2]C(=C[C]=C=O)CO'),
    E0 = (156.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,1185.73],'cm^-1')),
        HinderedRotor(inertia=(0.21575,'amu*angstrom^2'), symmetry=1, barrier=(4.96053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217623,'amu*angstrom^2'), symmetry=1, barrier=(5.00359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217423,'amu*angstrom^2'), symmetry=1, barrier=(4.99898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.60541,'amu*angstrom^2'), symmetry=1, barrier=(59.9035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823106,0.0758716,-0.000104834,8.76899e-08,-2.98037e-11,18882.5,30.0623], Tmin=(100,'K'), Tmax=(829.385,'K')), NASAPolynomial(coeffs=[7.09811,0.0362893,-1.6393e-05,3.05251e-09,-2.08042e-13,18162.2,2.89517], Tmin=(829.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]C(C=C=CO)=CO(29221)',
    structure = SMILES('[CH]C(C=C=CO)=CO'),
    E0 = (150.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3580,3650,1210,1345,900,1100,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.82674,'amu*angstrom^2'), symmetry=1, barrier=(42.0002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82983,'amu*angstrom^2'), symmetry=1, barrier=(42.0714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81446,'amu*angstrom^2'), symmetry=1, barrier=(41.718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82667,'amu*angstrom^2'), symmetry=1, barrier=(41.9988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.03928,0.108924,-0.00011773,5.78145e-08,-1.01977e-11,18335.1,35.5887], Tmin=(100,'K'), Tmax=(1701.95,'K')), NASAPolynomial(coeffs=[29.298,0.00119415,5.18231e-06,-1.2851e-09,9.2552e-14,11923.3,-124.115], Tmin=(1701.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C1[CH]C(=C[O])C1O(29222)',
    structure = SMILES('C=C1[CH]C(=C[O])C1O'),
    E0 = (52.6978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33545,0.0277006,9.27957e-05,-1.5705e-07,6.67799e-11,6462.57,26.6041], Tmin=(100,'K'), Tmax=(929.761,'K')), NASAPolynomial(coeffs=[25.1583,0.00449008,2.33811e-06,-4.78638e-10,2.01884e-14,-1394.02,-105.02], Tmin=(929.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.6978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = 'C=C1C=[C]C([O])C1O(29223)',
    structure = SMILES('C=C1C=[C]C([O])C1O'),
    E0 = (239.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06732,0.046874,1.33466e-05,-5.76295e-08,2.69112e-11,28891.9,26.4981], Tmin=(100,'K'), Tmax=(967.434,'K')), NASAPolynomial(coeffs=[18.8787,0.0160811,-5.34942e-06,1.0385e-09,-8.08977e-14,23440.3,-69.2058], Tmin=(967.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(3-Methylenecyclopentene) + radical(CC(C)OJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'OC=C1[CH][C]=COC1(29224)',
    structure = SMILES('OC=C1[CH][C]=COC1'),
    E0 = (57.0621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06687,0.0245395,0.00012907,-2.09186e-07,8.84329e-11,7005.61,19.5275], Tmin=(100,'K'), Tmax=(926.974,'K')), NASAPolynomial(coeffs=[31.1444,-0.00241429,6.28231e-06,-1.2037e-09,6.55619e-14,-2988.78,-147.129], Tmin=(926.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclohexane) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C=C(C=C=C=O)CO(29225)',
    structure = SMILES('C=C(C=C=C=O)CO'),
    E0 = (20.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28551,0.0525272,-3.8363e-05,9.11606e-09,1.18281e-12,2574.4,11.4247], Tmin=(100,'K'), Tmax=(1016.24,'K')), NASAPolynomial(coeffs=[13.7311,0.017109,-6.11226e-06,1.09751e-09,-7.66957e-14,-655.773,-52.2657], Tmin=(1016.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = '[CH]C(=CO)C[C]=C[O](25677)',
    structure = SMILES('[CH]C(=CO)C[C]=C[O]'),
    E0 = (382.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.27202,0.0767641,-3.5312e-05,-2.89789e-08,2.26202e-11,46125.6,31.626], Tmin=(100,'K'), Tmax=(917.058,'K')), NASAPolynomial(coeffs=[25.3738,0.010315,-9.02425e-07,4.64301e-12,-1.44907e-15,39512.3,-100.291], Tmin=(917.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]C([CH]C=C[O])=CO(25674)',
    structure = SMILES('[CH]C([CH]C=C[O])=CO'),
    E0 = (245.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.07515,0.0598143,2.70559e-05,-1.03071e-07,5.1553e-11,29750.9,29.5776], Tmin=(100,'K'), Tmax=(908.965,'K')), NASAPolynomial(coeffs=[28.8136,0.00363593,3.77173e-06,-9.21448e-10,5.99719e-14,21622.8,-122.298], Tmin=(908.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(C=[C]C[O])=CO(29226)',
    structure = SMILES('[CH]C(C=[C]C[O])=CO'),
    E0 = (492.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.164725,0.0857769,-8.9748e-05,4.87036e-08,-1.04252e-11,59402,31.4067], Tmin=(100,'K'), Tmax=(1142.57,'K')), NASAPolynomial(coeffs=[17.4323,0.0241719,-8.87119e-06,1.51374e-09,-9.98273e-14,55380.8,-55.8342], Tmin=(1142.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(492.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'CC([CH][C]=[C][O])=CO(29227)',
    structure = SMILES('CC([CH][C]=[C][O])=CO'),
    E0 = (352.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,308.609,308.61,308.611],'cm^-1')),
        HinderedRotor(inertia=(0.315365,'amu*angstrom^2'), symmetry=1, barrier=(21.314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.315365,'amu*angstrom^2'), symmetry=1, barrier=(21.314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.315366,'amu*angstrom^2'), symmetry=1, barrier=(21.314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51955,'amu*angstrom^2'), symmetry=1, barrier=(102.695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.254205,0.066243,-2.51245e-05,-3.3993e-08,2.37521e-11,42592,31.4276], Tmin=(100,'K'), Tmax=(909.215,'K')), NASAPolynomial(coeffs=[23.587,0.00709252,6.95634e-07,-3.04863e-10,2.05025e-14,36551.1,-88.8067], Tmin=(909.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C(O)C=C=C[O](25844)',
    structure = SMILES('C=[C]C(O)C=C=C[O]'),
    E0 = (222.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.889253,'amu*angstrom^2'), symmetry=1, barrier=(20.4457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888425,'amu*angstrom^2'), symmetry=1, barrier=(20.4266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888295,'amu*angstrom^2'), symmetry=1, barrier=(20.4236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4502.04,'J/mol'), sigma=(7.09037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=703.21 K, Pc=28.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0160647,0.0765635,-7.56417e-05,3.67375e-08,-6.91216e-12,26854.1,32.81], Tmin=(100,'K'), Tmax=(1306.09,'K')), NASAPolynomial(coeffs=[19.9411,0.0155418,-5.5605e-06,9.66045e-10,-6.51257e-14,21649.3,-68.6373], Tmin=(1306.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ)"""),
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
    label = 'C=[C]C=C=C[O](23378)',
    structure = SMILES('C=[C]C=C=C[O]'),
    E0 = (366.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63597,'amu*angstrom^2'), symmetry=1, barrier=(37.6142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442946,0.0581895,-5.88988e-05,2.83071e-08,-4.99578e-12,44184.9,22.5047], Tmin=(100,'K'), Tmax=(1656.64,'K')), NASAPolynomial(coeffs=[16.9916,0.00540583,5.07716e-07,-2.72871e-10,2.24426e-14,40461.9,-60.3744], Tmin=(1656.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C1([CH]O)C=C1C=O(29107)',
    structure = SMILES('[CH2]C1([CH]O)C=C1C=O'),
    E0 = (302.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.62604,0.0805137,-9.38063e-05,5.83392e-08,-1.49379e-11,36492.7,26.2699], Tmin=(100,'K'), Tmax=(935.104,'K')), NASAPolynomial(coeffs=[11.8606,0.0324568,-1.67182e-05,3.38059e-09,-2.44718e-13,34391.6,-27.1766], Tmin=(935.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopropene) + radical(CCsJOH) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C(C#CC=O)=CO(29228)',
    structure = SMILES('[CH2]C(C#CC=O)=CO'),
    E0 = (42.9971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2100,2250,500,550,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.21406,'amu*angstrom^2'), symmetry=1, barrier=(27.9136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21496,'amu*angstrom^2'), symmetry=1, barrier=(27.9344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21526,'amu*angstrom^2'), symmetry=1, barrier=(27.9411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21573,'amu*angstrom^2'), symmetry=1, barrier=(27.9521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.367044,0.0669906,-4.21934e-05,-7.572e-09,1.17331e-11,5314,26.4691], Tmin=(100,'K'), Tmax=(939.988,'K')), NASAPolynomial(coeffs=[21.3399,0.00982173,-2.15533e-06,3.37188e-10,-2.62454e-14,-46.0475,-80.9533], Tmin=(939.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.9971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=CC=O)=CO(29229)',
    structure = SMILES('[CH2]C([C]=CC=O)=CO'),
    E0 = (74.5446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.819289,0.0901333,-9.80427e-05,5.0688e-08,-9.89506e-12,9152.7,28.7117], Tmin=(100,'K'), Tmax=(1369.4,'K')), NASAPolynomial(coeffs=[24.5288,0.00922133,-1.88841e-06,2.13283e-10,-1.14464e-14,2854.54,-99.1951], Tmin=(1369.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.5446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([CH]O)=CC=C=O(29230)',
    structure = SMILES('[CH2]C([CH]O)=CC=C=O'),
    E0 = (35.5485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736179,0.0639467,-5.22918e-05,2.21304e-08,-3.76559e-12,4399.49,30.5043], Tmin=(100,'K'), Tmax=(1402.55,'K')), NASAPolynomial(coeffs=[14.8917,0.0235763,-9.11709e-06,1.60869e-09,-1.07707e-13,428.664,-42.5769], Tmin=(1402.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.5485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJO) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(=[C]O)C=CC=O(29231)',
    structure = SMILES('[CH2]C(=[C]O)C=CC=O'),
    E0 = (115.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.969793,'amu*angstrom^2'), symmetry=1, barrier=(22.2975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970443,'amu*angstrom^2'), symmetry=1, barrier=(22.3124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970049,'amu*angstrom^2'), symmetry=1, barrier=(22.3033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970591,'amu*angstrom^2'), symmetry=1, barrier=(22.3158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.283835,0.0823935,-8.52949e-05,4.25525e-08,-8.16029e-12,14030.8,30.3239], Tmin=(100,'K'), Tmax=(1287,'K')), NASAPolynomial(coeffs=[21.9467,0.0133006,-4.76685e-06,8.38827e-10,-5.73779e-14,8308.69,-82.5348], Tmin=(1287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJO) + radical(Allyl_P)"""),
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
    label = 'O=CC1=CC(=CO)C1(29089)',
    structure = SMILES('O=CC1=CC(=CO)C1'),
    E0 = (-119.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510143,0.0538963,1.51677e-05,-7.5177e-08,3.69685e-11,-14271,23.523], Tmin=(100,'K'), Tmax=(943.966,'K')), NASAPolynomial(coeffs=[25.3334,0.00611362,-1.21453e-07,4.23869e-11,-1.40203e-14,-21515,-108.35], Tmin=(943.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene)"""),
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
    label = '[C]=CC([CH2])=CO(28232)',
    structure = SMILES('[C]=CC([CH2])=CO'),
    E0 = (557.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.22113,'amu*angstrom^2'), symmetry=1, barrier=(28.0763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22167,'amu*angstrom^2'), symmetry=1, barrier=(28.0887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21961,'amu*angstrom^2'), symmetry=1, barrier=(28.0413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.541262,0.0716231,-7.75329e-05,3.81019e-08,-6.74676e-12,67196.3,25.4375], Tmin=(100,'K'), Tmax=(1681.01,'K')), NASAPolynomial(coeffs=[21.2348,0.00044384,3.25953e-06,-7.91454e-10,5.65178e-14,62610.9,-82.7926], Tmin=(1681.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]C(C=CC=O)=CO(29232)',
    structure = SMILES('[CH]C(C=CC=O)=CO'),
    E0 = (94.7344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.87949,'amu*angstrom^2'), symmetry=1, barrier=(43.2132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87826,'amu*angstrom^2'), symmetry=1, barrier=(43.1848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87812,'amu*angstrom^2'), symmetry=1, barrier=(43.1816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87949,'amu*angstrom^2'), symmetry=1, barrier=(43.2132,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254634,0.0803617,-5.77914e-05,6.14456e-09,6.28698e-12,11558.9,27.9115], Tmin=(100,'K'), Tmax=(980.357,'K')), NASAPolynomial(coeffs=[22.0806,0.0180176,-6.44741e-06,1.18167e-09,-8.55371e-14,5796.24,-86.4556], Tmin=(980.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.7344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
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
    label = 'C=C1C=C(C=O)C1O(29105)',
    structure = SMILES('C=C1C=C(C=O)C1O'),
    E0 = (-88.9148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.837827,0.055186,-1.4627e-05,-2.40306e-08,1.33774e-11,-10567.6,24.1821], Tmin=(100,'K'), Tmax=(1027.24,'K')), NASAPolynomial(coeffs=[18.061,0.0195255,-8.41346e-06,1.69861e-09,-1.27451e-13,-15763,-67.4377], Tmin=(1027.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.9148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    E0 = (72.2858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (220.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (436.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (331.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (298.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (317.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (296.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (327.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (231.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (449.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (137.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (347.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (144.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (528.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (425.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (456.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (508.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (524.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (467.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (176.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (259.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (202.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (204.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (203.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (218.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (97.2591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (368.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (511.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (255.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (294.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (466.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (274.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (322.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (451.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (521.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (193.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (380.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (439.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (361.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (216.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (573.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (688.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (503.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (256.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (285.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (314.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (320.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (368.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (220.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (301.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (197.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (239.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (218.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (80.6538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (97.2591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (460.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (270.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (517.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (377.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (316.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (606.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (303.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (268.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (287.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (259.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (235.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (518.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (105.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (80.5701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (624.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (498.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (111.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (80.1936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['C=C=CO(12571)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[O]C=[C]C1CC1=CO(29139)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.32e+10,'s^-1'), n=0.35, Ea=(148.079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 147.0 to 148.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C1([CH]O)C=C=CO1(29192)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(364.525,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMM;doublebond_intra;radadd_intra] for rate rule [R6_SMM;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2]C(C=C=C=O)=CO(29193)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C=C[O](8556)', 'C=C=CO(12571)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.04,'cm^3/(mol*s)'), n=3.05, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2838 used for Ca_Cds-HH;CdsJ=Cdd
Exact match found for rate rule [Ca_Cds-HH;CdsJ=Cdd]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C([CH]O)=CC#CO(29194)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC(=C=C=C[O])[CH]O(29195)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_double;Cs_H_out_2H] for rate rule [R3H_SS_2Cd;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['CC(=[C]O)C=C=C[O](29196)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['CC(C=O)=C[C]=C[O](25373)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(493835,'s^-1'), n=1.76395, Ea=(159.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_2H;XH_out] for rate rule [R4H_SDS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=C=C=CO)[CH]O(29197)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.63625e+10,'s^-1'), n=1.0925, Ea=(294.328,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R4H_MMS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC([CH][C]=C=O)=CO(29198)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(=[C]O)C=C=CO(29199)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;XH_out] for rate rule [R6H;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=C=CO)=C[O](25379)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(627096,'s^-1'), n=1.03067, Ea=(72.2138,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;O_rad_out;XH_out] for rate rule [R7H;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OH(D)(132)', '[CH]C(=C)C=C=C[O](22618)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH2]C(C=C=C[O])=C[O](25381)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C=C[O](8556)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH2]C(=C=C=C[O])[CH]O(29200)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH2]C(=[C]O)C=C=C[O](29201)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2]C([CH][C]=C=O)=CO(29202)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[O]C=C=C[C]1CC1O(29203)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.03537e+13,'s^-1'), n=-0.296394, Ea=(104.177,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C(=CO)C=C1[CH]O1(29204)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[O]C=C1[CH]C(=CO)C1(29205)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C(=CO)C1[C]=CO1(29206)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(131.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 129.7 to 132.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[O]C1[C]=CC(=CO)C1(29207)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.4946e+10,'s^-1'), n=0.314866, Ea=(130.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 127.0 to 130.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C1=C[C]=COC1O(29208)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0.19, Ea=(146.44,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMM;doublebond_intra;radadd_intra] for rate rule [R6_SMM;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['CC(C=C=C=O)=CO(29209)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(=C[C]=C[O])C[O](25394)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([C]=C=C[O])[CH]O(25400)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(=[C]C=C[O])[CH]O(29210)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([CH]C=[C][O])=CO(29211)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(=[C]O)C[C]=C[O](29212)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(=[C]O)[CH]C=C[O](29213)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(=C[O])C[C]=C[O](25393)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH]O)[CH][C]=C=O(25405)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=[C]O)C=[C]C[O](29214)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([CH]C=C[O])=C[O](25399)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(=C[C]=[C][O])CO(29215)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(C=[C]C[O])=C[O](25402)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C(C=O)C=C=C[O](22596)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(T)(28)', '[O]C=C=C[C]=CO(29216)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O(T)(63)', 'C#CC=C([CH2])[CH]O(27790)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(8)', '[CH]C(C=C=C[O])=CO(29217)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['C=C1C(O)C1[C]=C[O](29165)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(8.15385e+14,'s^-1'), n=-0.537569, Ea=(184.065,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 182.4 to 184.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction45',
    reactants = ['H(8)', 'C=C(C=O)C=C=C[O](25371)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2834 used for Od_CO-CdH;HJ
Exact match found for rate rule [Od_CO-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C(C=C=C[O])C[O](26530)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(=C=C=C[O])CO(29218)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_double;Cs_H_out_1H] for rate rule [R3H_SS_2Cd;Cd_rad_out_double;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C(C=C=C[O])CO(29219)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C(=C[C]=C=O)CO(29220)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(31338.3,'s^-1'), n=2.07906, Ea=(64.7147,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5H;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]C(C=C=CO)=CO(29221)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['C=C1[CH]C(=C[O])C1O(29222)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(3.6836e+08,'s^-1'), n=0.948854, Ea=(125.32,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['C=C1C=[C]C([O])C1O(29223)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(166.925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra_csHNd] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_csHO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 162.5 to 166.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['OC=C1[CH][C]=COC1(29224)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0.19, Ea=(146.44,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMM;doublebond_intra;radadd_intra] for rate rule [R6_SMM;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['C=C(C=C=C=O)CO(29225)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['C=C(C=O)C=C=CO(25391)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]C(=CO)C[C]=C[O](25677)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH]C([CH]C=C[O])=CO(25674)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH]C(C=[C]C[O])=CO(29226)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction59',
    reactants = ['CC([CH][C]=[C][O])=CO(29227)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH]O(5471)', 'C=[C]C=C=C[O](23378)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C1([CH]O)C=C1C=O(29107)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1.126e+14,'s^-1'), n=-0.355, Ea=(230.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingleDe] for rate rule [R4_D_D;doublebond_intra_HNd;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['H(8)', '[CH2]C(C#CC=O)=CO(29228)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2][C]=CO(18753)', 'C#CC=O(21959)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C([C]=CC=O)=CO(29229)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(7.7652e+07,'s^-1'), n=1.65613, Ea=(187.464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_D;Cd_rad_out_single;Cd_H_out_singleDe] + [R2H_D;Cd_rad_out_singleDe;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C([CH]O)=CC=C=O(29230)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH2]C(=[C]O)C=CC=O(29231)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['[CH2]C(C=O)=CC=C[O](25414)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;XH_out] for rate rule [R5H_DSMS;Cd_rad_out_singleDe;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['O=CC1=CC(=CO)C1(29089)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH]=O(373)', '[C]=CC([CH2])=CO(28232)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]C(C=CC=O)=CO(29232)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['C=C(C=O)C=CC=O(25418)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH2]C(C=C=C[O])=CO(25407)'],
    products = ['C=C1C=C(C=O)C1O(29105)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeO;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

network(
    label = '5259',
    isomers = [
        '[CH2]C(C=C=C[O])=CO(25407)',
    ],
    reactants = [
        ('C=C=CO(12571)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5259',
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

