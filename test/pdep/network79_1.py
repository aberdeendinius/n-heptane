species(
    label = '[CH]=CC[CH]C(151)',
    structure = SMILES('[CH]=CC[CH]C'),
    E0 = (401.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.00455381,'amu*angstrom^2'), symmetry=1, barrier=(5.17142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225503,'amu*angstrom^2'), symmetry=1, barrier=(5.18475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226422,'amu*angstrom^2'), symmetry=1, barrier=(5.20589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05382,0.0394038,-2.00486e-05,4.87108e-09,-4.71722e-13,48331,22.4667], Tmin=(100,'K'), Tmax=(2304.45,'K')), NASAPolynomial(coeffs=[13.2679,0.0199387,-7.37854e-06,1.2057e-09,-7.40814e-14,43162.5,-40.9971], Tmin=(2304.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJC)"""),
)

species(
    label = 'C#C(139)',
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
    label = 'C=CC(37)',
    structure = SMILES('C=CC'),
    E0 = (6.12372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.597443,'amu*angstrom^2'), symmetry=1, barrier=(13.7364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30977,0.00827491,3.37717e-05,-4.3931e-08,1.58773e-11,767.476,9.64349], Tmin=(100,'K'), Tmax=(988,'K')), NASAPolynomial(coeffs=[5.41204,0.0172866,-6.51359e-06,1.20323e-09,-8.55924e-14,-503.177,-4.80153], Tmin=(988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.12372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'H(7)',
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
    label = '[CH]=CC=CC(181)',
    structure = SMILES('[CH]=CC=CC'),
    E0 = (304.992,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680],'cm^-1')),
        HinderedRotor(inertia=(0.867878,'amu*angstrom^2'), symmetry=1, barrier=(19.9542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87242,'amu*angstrom^2'), symmetry=1, barrier=(20.0587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93303,0.0363756,-1.10322e-06,-2.51399e-08,1.24819e-11,36764.6,17.6229], Tmin=(100,'K'), Tmax=(973.629,'K')), NASAPolynomial(coeffs=[11.5283,0.0184065,-6.46811e-06,1.16265e-09,-8.23224e-14,33879.4,-33.6335], Tmin=(973.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCC=C(182)',
    structure = SMILES('[CH]=CCC=C'),
    E0 = (335.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,249.876],'cm^-1')),
        HinderedRotor(inertia=(0.328692,'amu*angstrom^2'), symmetry=1, barrier=(14.5402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327888,'amu*angstrom^2'), symmetry=1, barrier=(14.5399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15389,0.0316791,6.94394e-06,-2.91267e-08,1.26951e-11,40467.3,19.8536], Tmin=(100,'K'), Tmax=(1001.64,'K')), NASAPolynomial(coeffs=[10.0557,0.0208756,-7.95414e-06,1.47296e-09,-1.04739e-13,37843.3,-23.4776], Tmin=(1001.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C#CC[CH]C(183)',
    structure = SMILES('C#CC[CH]C'),
    E0 = (320.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.196151,'amu*angstrom^2'), symmetry=1, barrier=(4.50989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195987,'amu*angstrom^2'), symmetry=1, barrier=(4.50612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0435023,'amu*angstrom^2'), symmetry=1, barrier=(64.2953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23907,0.0380811,-2.41036e-05,9.47658e-09,-1.71378e-12,38580.8,19.2199], Tmin=(100,'K'), Tmax=(1210.43,'K')), NASAPolynomial(coeffs=[5.6734,0.0267321,-1.00397e-05,1.73069e-09,-1.13976e-13,37749.3,1.9953], Tmin=(1210.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJC)"""),
)

species(
    label = '[CH]=[CH](140)',
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
    label = '[CH2][CH]C(39)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (279.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00418542,'amu*angstrom^2'), symmetry=1, barrier=(6.91845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00418543,'amu*angstrom^2'), symmetry=1, barrier=(6.91844,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25506,0.0137284,1.00539e-05,-1.43791e-08,4.38762e-12,33590.4,14.1736], Tmin=(100,'K'), Tmax=(1201.85,'K')), NASAPolynomial(coeffs=[3.74305,0.0203098,-8.40112e-06,1.53862e-09,-1.05138e-13,32880.5,9.26418], Tmin=(1201.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH]C=CCC(177)',
    structure = SMILES('[CH]C=CCC'),
    E0 = (318.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,383.524,384.989,385.082,387.907],'cm^-1')),
        HinderedRotor(inertia=(0.492463,'amu*angstrom^2'), symmetry=1, barrier=(51.4859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.532443,'amu*angstrom^2'), symmetry=1, barrier=(51.5194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48962,'amu*angstrom^2'), symmetry=1, barrier=(51.4438,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95737,0.0362808,9.16369e-06,-2.9438e-08,1.16674e-11,38344,20.3268], Tmin=(100,'K'), Tmax=(1046.98,'K')), NASAPolynomial(coeffs=[8.10496,0.0322442,-1.29197e-05,2.36788e-09,-1.64298e-13,35990.7,-14.705], Tmin=(1046.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CCC[CH2](184)',
    structure = SMILES('[CH]=CCC[CH2]'),
    E0 = (412.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.137724,'amu*angstrom^2'), symmetry=1, barrier=(3.16656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142612,'amu*angstrom^2'), symmetry=1, barrier=(3.27892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0401954,'amu*angstrom^2'), symmetry=1, barrier=(14.3342,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84058,0.0407678,-1.49453e-05,-5.21373e-09,3.75887e-12,49640.4,22.2665], Tmin=(100,'K'), Tmax=(1101.28,'K')), NASAPolynomial(coeffs=[9.60315,0.0248089,-9.87415e-06,1.80514e-09,-1.24698e-13,47188.7,-19.3011], Tmin=(1101.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C[CH]C(185)',
    structure = SMILES('C=[C]C[CH]C'),
    E0 = (391.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.41566,0.0369675,-1.69643e-05,3.43217e-09,-2.57685e-13,47199.8,21.2179], Tmin=(100,'K'), Tmax=(2427.98,'K')), NASAPolynomial(coeffs=[16.5625,0.0166401,-6.24654e-06,9.9465e-10,-5.87348e-14,39452.1,-61.39], Tmin=(2427.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CCC(186)',
    structure = SMILES('[CH]=[C]CCC'),
    E0 = (444.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,1685,370,352.047],'cm^-1')),
        HinderedRotor(inertia=(0.102332,'amu*angstrom^2'), symmetry=1, barrier=(9.00278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101878,'amu*angstrom^2'), symmetry=1, barrier=(9.00673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101757,'amu*angstrom^2'), symmetry=1, barrier=(9.01071,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80017,0.0442028,-2.73482e-05,8.62695e-09,-1.12131e-12,53559.8,21.1032], Tmin=(100,'K'), Tmax=(1729.9,'K')), NASAPolynomial(coeffs=[10.9041,0.0231521,-9.09516e-06,1.59261e-09,-1.04726e-13,50410,-27.8072], Tmin=(1729.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C[CH]C(187)',
    structure = SMILES('[CH2]C=C[CH]C'),
    E0 = (240.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18177,0.028357,2.70944e-05,-5.14677e-08,2.0569e-11,28948.9,17.5848], Tmin=(100,'K'), Tmax=(990.215,'K')), NASAPolynomial(coeffs=[10.237,0.0240424,-9.12509e-06,1.70242e-09,-1.22293e-13,25969.8,-28.1847], Tmin=(990.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC=C(188)',
    structure = SMILES('[CH2][CH]CC=C'),
    E0 = (359.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13435,0.0369817,-1.51233e-05,1.33711e-09,3.7705e-13,43295.2,23.5599], Tmin=(100,'K'), Tmax=(1532.35,'K')), NASAPolynomial(coeffs=[8.68443,0.0259362,-1.02359e-05,1.7885e-09,-1.17141e-13,40577.2,-13.1548], Tmin=(1532.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH]C=C[CH]C(189)',
    structure = SMILES('[CH]C=C[CH]C'),
    E0 = (459.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,484.313,485.094,485.318,486.16],'cm^-1')),
        HinderedRotor(inertia=(0.309021,'amu*angstrom^2'), symmetry=1, barrier=(51.6012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308971,'amu*angstrom^2'), symmetry=1, barrier=(51.5978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308845,'amu*angstrom^2'), symmetry=1, barrier=(51.5774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1723,0.0308697,1.87068e-05,-3.85804e-08,1.48646e-11,55308.9,18.4482], Tmin=(100,'K'), Tmax=(1025.92,'K')), NASAPolynomial(coeffs=[8.03499,0.029944,-1.20077e-05,2.21714e-09,-1.55209e-13,52951.7,-15.6112], Tmin=(1025.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC[CH][CH2](190)',
    structure = SMILES('[CH]=CC[CH][CH2]'),
    E0 = (606.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1500.53],'cm^-1')),
        HinderedRotor(inertia=(0.0958967,'amu*angstrom^2'), symmetry=1, barrier=(2.20485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965059,'amu*angstrom^2'), symmetry=1, barrier=(2.21886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097205,'amu*angstrom^2'), symmetry=1, barrier=(2.23494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43023,0.0365998,-1.99417e-05,5.48319e-09,-6.43181e-13,72997.5,22.8669], Tmin=(100,'K'), Tmax=(1732.14,'K')), NASAPolynomial(coeffs=[7.07421,0.0258755,-1.06546e-05,1.90871e-09,-1.27271e-13,71388.7,-2.08879], Tmin=(1732.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P) + radical(RCCJC)"""),
)

species(
    label = '[CH]=[C]C[CH]C(191)',
    structure = SMILES('[CH]=[C]C[CH]C'),
    E0 = (639.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,3025,407.5,1350,352.5,232.393],'cm^-1')),
        HinderedRotor(inertia=(0.154267,'amu*angstrom^2'), symmetry=1, barrier=(5.8792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152567,'amu*angstrom^2'), symmetry=1, barrier=(5.88518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00275597,'amu*angstrom^2'), symmetry=1, barrier=(5.88659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91207,0.0302162,2.78285e-05,-1.26167e-07,1.1567e-10,76896.3,20.027], Tmin=(100,'K'), Tmax=(409.625,'K')), NASAPolynomial(coeffs=[3.84828,0.0313221,-1.3748e-05,2.57466e-09,-1.7801e-13,76733.6,15.2964], Tmin=(409.625,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=CC(192)',
    structure = SMILES('C=CC=CC'),
    E0 = (57.8956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00727,0.0328459,1.55855e-05,-4.25745e-08,1.84259e-11,7044.82,16.9534], Tmin=(100,'K'), Tmax=(972.32,'K')), NASAPolynomial(coeffs=[11.2869,0.0212416,-7.50361e-06,1.3618e-09,-9.72233e-14,3984.25,-34.0139], Tmin=(972.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.8956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCC=C(193)',
    structure = SMILES('C=CCC=C'),
    E0 = (88.7498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23067,0.0281207,2.37269e-05,-4.66715e-08,1.86798e-11,10747.4,18.4818], Tmin=(100,'K'), Tmax=(993.874,'K')), NASAPolynomial(coeffs=[9.79588,0.0237415,-9.00715e-06,1.6762e-09,-1.19976e-13,7956.11,-24.4466], Tmin=(993.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.7498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]=CC([CH2])C(152)',
    structure = SMILES('[CH]=CC([CH2])C'),
    E0 = (403.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(3.49271,'amu*angstrom^2'), symmetry=1, barrier=(80.3042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693181,'amu*angstrom^2'), symmetry=1, barrier=(15.9376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00226235,'amu*angstrom^2'), symmetry=1, barrier=(15.9621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3076.04,'J/mol'), sigma=(5.56648,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=480.47 K, Pc=40.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91876,0.0379709,-2.37229e-06,-2.30889e-08,1.17463e-11,48660.5,21.7321], Tmin=(100,'K'), Tmax=(953.434,'K')), NASAPolynomial(coeffs=[10.1215,0.0224181,-7.57623e-06,1.29768e-09,-8.83996e-14,46239.1,-21.9457], Tmin=(953.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = 'CC1C=CC1(178)',
    structure = SMILES('CC1C=CC1'),
    E0 = (112.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71863,0.00906983,8.67845e-05,-1.17242e-07,4.49737e-11,13620,15.2056], Tmin=(100,'K'), Tmax=(953.248,'K')), NASAPolynomial(coeffs=[11.8044,0.0195073,-6.05665e-06,1.13153e-09,-8.70693e-14,9681.42,-39.766], Tmin=(953.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]C(29)',
    structure = SMILES('[CH]C'),
    E0 = (351.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,431.535,1804.51],'cm^-1')),
        HinderedRotor(inertia=(0.000906356,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73285,-0.000244773,3.59198e-05,-4.44289e-08,1.65887e-11,42287.5,7.078], Tmin=(100,'K'), Tmax=(940.479,'K')), NASAPolynomial(coeffs=[5.42969,0.0081677,-2.4253e-06,4.22642e-10,-3.09417e-14,41277.1,-4.6789], Tmin=(940.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C=C(163)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,192.655,193.544,193.915],'cm^-1')),
        HinderedRotor(inertia=(1.88068,'amu*angstrom^2'), symmetry=1, barrier=(50.3487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32096,0.00806329,3.46645e-05,-4.52343e-08,1.64854e-11,45350.1,10.7121], Tmin=(100,'K'), Tmax=(975.253,'K')), NASAPolynomial(coeffs=[5.21066,0.0176207,-6.65616e-06,1.20944e-09,-8.49962e-14,44158.4,-2.57721], Tmin=(975.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]CC=[CH](164)',
    structure = SMILES('[CH]CC=[CH]'),
    E0 = (678.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,673.753,674.423,675.591],'cm^-1')),
        HinderedRotor(inertia=(0.153602,'amu*angstrom^2'), symmetry=1, barrier=(3.53162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152901,'amu*angstrom^2'), symmetry=1, barrier=(3.51549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51662,0.0269418,-5.60112e-06,-1.14302e-08,6.17492e-12,81697.9,16.8097], Tmin=(100,'K'), Tmax=(1000.99,'K')), NASAPolynomial(coeffs=[8.88798,0.0141078,-5.28978e-06,9.63822e-10,-6.77525e-14,79789.8,-17.0944], Tmin=(1000.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(678.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH3](10)',
    structure = SMILES('[CH3]'),
    E0 = (135.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([570.572,1408.13,1408.49,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91547,0.00184154,3.48742e-06,-3.32748e-09,8.49957e-13,16285.6,0.351741], Tmin=(100,'K'), Tmax=(1337.63,'K')), NASAPolynomial(coeffs=[3.54146,0.00476787,-1.82148e-06,3.28877e-10,-2.22546e-14,16224,1.66035], Tmin=(1337.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""),
)

species(
    label = '[CH]=CC[C]C(194)',
    structure = SMILES('[CH]=CC[C]C'),
    E0 = (655.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,180,759.276],'cm^-1')),
        HinderedRotor(inertia=(0.844404,'amu*angstrom^2'), symmetry=1, barrier=(19.4145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386829,'amu*angstrom^2'), symmetry=1, barrier=(8.89396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.847008,'amu*angstrom^2'), symmetry=1, barrier=(19.4744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74879,0.0431294,-2.4222e-05,2.60598e-09,1.59284e-12,78865.9,20.714], Tmin=(100,'K'), Tmax=(1098.69,'K')), NASAPolynomial(coeffs=[10.6138,0.0211744,-8.33718e-06,1.51654e-09,-1.04586e-13,76295,-25.7235], Tmin=(1098.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=CC[CH]C(195)',
    structure = SMILES('[C]=CC[CH]C'),
    E0 = (712.246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,217.66,1650.94],'cm^-1')),
        HinderedRotor(inertia=(0.141119,'amu*angstrom^2'), symmetry=1, barrier=(5.10957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210768,'amu*angstrom^2'), symmetry=1, barrier=(5.0633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318046,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86342,0.0342106,-1.63058e-05,3.43633e-09,-2.79402e-13,85693.8,19.4493], Tmin=(100,'K'), Tmax=(2650.63,'K')), NASAPolynomial(coeffs=[12.206,0.0201119,-8.32725e-06,1.42961e-09,-9.01341e-14,80741.1,-34.7304], Tmin=(2650.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(712.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(RCCJC)"""),
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
    E0 = (401.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (524.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (550.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (552.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (554.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (521.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (508.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (563.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (506.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (586.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (593.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (434.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (815.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (671.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (818.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (850.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (464.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (426.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (563.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (409.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (762.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (848.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (866.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (924.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC[CH]C(151)'],
    products = ['C#C(139)', 'C=CC(37)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(7)', '[CH]=CC=CC(181)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(105.529,'m^3/(mol*s)'), n=1.6629, Ea=(8.08712,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(7)', '[CH]=CCC=C(182)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(7)', 'C#CC[CH]C(183)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.276e+10,'cm^3/(mol*s)'), n=1.103, Ea=(20.5476,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 139 used for Ct-Cs_Ct-H;HJ
Exact match found for rate rule [Ct-Cs_Ct-H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[CH](140)', 'C=CC(37)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0124089,'m^3/(mol*s)'), n=2.46568, Ea=(12.4666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C#C(139)', '[CH2][CH]C(39)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(46.4627,'m^3/(mol*s)'), n=1.51997, Ea=(27.4714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C=CCC(177)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.9e+10,'s^-1'), n=0.75, Ea=(190.79,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 159 used for R2H_S;C_rad_out_H/Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=CCC[CH2](184)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC[CH]C(151)'],
    products = ['C=[C]C[CH]C(185)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CCC(186)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC[CH]C(151)'],
    products = ['[CH2]C=C[CH]C(187)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC[CH]C(151)'],
    products = ['[CH2][CH]CC=C(188)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[CH](140)', '[CH2][CH]C(39)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(7)', '[CH]C=C[CH]C(189)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(7)', '[CH]=CC[CH][CH2](190)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(7)', '[CH]=[C]C[CH]C(191)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC[CH]C(151)'],
    products = ['C=CC=CC(192)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC[CH]C(151)'],
    products = ['C=CCC=C(193)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC([CH2])C(152)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC[CH]C(151)'],
    products = ['CC1C=CC1(178)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(29)', '[CH]C=C(163)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]CC=[CH](164)', '[CH3](10)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(7)', '[CH]=CC[C]C(194)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[C]=CC[CH]C(195)', 'H(7)'],
    products = ['[CH]=CC[CH]C(151)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '79',
    isomers = [
        '[CH]=CC[CH]C(151)',
    ],
    reactants = [
        ('C#C(139)', 'C=CC(37)'),
    ],
    bathGas = {
        'Ne': 0.25,
        'He': 0.25,
        'N2': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '79',
    Tmin = (300,'K'),
    Tmax = (3000,'K'),
    Tcount = 8,
    Tlist = ([302.617,324.619,374.997,470.374,649.057,1000.02,1706.11,2761.25],'K'),
    Pmin = (0.001,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.00132544,0.0107284,0.316228,9.32101,75.4469],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

