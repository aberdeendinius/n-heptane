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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.694207,0.0803134,-8.09633e-05,3.89328e-08,-6.96019e-12,57599.4,31.4631], Tmin=(100,'K'), Tmax=(1590.54,'K')), NASAPolynomial(coeffs=[22.5117,0.00849225,-5.35427e-07,-9.9197e-11,1.1194e-14,51920.1,-85.9085], Tmin=(1590.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P) + radical(Cds_S)"""),
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
    label = '[CH]=C1CC1([CH2])[CH]O(28742)',
    structure = SMILES('[CH]=C1CC1([CH2])[CH]O'),
    E0 = (595.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.349984,0.0715004,-6.81919e-05,3.26422e-08,-6.12197e-12,71706.9,27.138], Tmin=(100,'K'), Tmax=(1300.93,'K')), NASAPolynomial(coeffs=[17.6498,0.0183082,-6.86002e-06,1.21234e-09,-8.20795e-14,67205.7,-60.8747], Tmin=(1300.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Neopentyl) + radical(Cds_P) + radical(CCsJOH)"""),
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
    label = '[CH]C=CC([CH2])=CO(28744)',
    structure = SMILES('[CH]C=CC([CH2])=CO'),
    E0 = (333.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.161647,0.0655937,-8.8067e-06,-5.31581e-08,3.05363e-11,40290.1,25.4285], Tmin=(100,'K'), Tmax=(916.865,'K')), NASAPolynomial(coeffs=[23.9154,0.0114095,-1.05539e-06,2.55007e-11,-3.46247e-15,33856,-98.4417], Tmin=(916.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#C[CH][C](C)[CH]O(28745)',
    structure = SMILES('C#C[CH][C](C)[CH]O'),
    E0 = (412.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.780936,0.0714448,-7.63767e-05,4.11304e-08,-7.07113e-12,49695.9,25.9824], Tmin=(100,'K'), Tmax=(783.615,'K')), NASAPolynomial(coeffs=[11.9881,0.0251803,-8.76432e-06,1.42987e-09,-9.086e-14,47603.4,-27.497], Tmin=(783.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJ(C)CO) + radical(Sec_Propargyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=[C]CC(C)=[C]O(28746)',
    structure = SMILES('[CH]=[C]CC(C)=[C]O'),
    E0 = (565.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.473494,0.0735421,-7.89856e-05,4.39732e-08,-9.63252e-12,68154.5,30.1023], Tmin=(100,'K'), Tmax=(1118.08,'K')), NASAPolynomial(coeffs=[15.3951,0.02016,-7.36975e-06,1.27206e-09,-8.47528e-14,64817.8,-43.5512], Tmin=(1118.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S) + radical(Cds_P)"""),
)

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
    label = '[CH]=CCC([CH2])=[C]O(28747)',
    structure = SMILES('[CH]=CCC([CH2])=[C]O'),
    E0 = (479.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.505461,0.0664512,-4.56222e-05,1.26089e-09,7.40452e-12,57774.8,29.7947], Tmin=(100,'K'), Tmax=(942.586,'K')), NASAPolynomial(coeffs=[18.6198,0.0150403,-4.32506e-06,7.08693e-10,-4.9462e-14,53228.9,-62.5252], Tmin=(942.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]CC(C)=C[O](19967)',
    structure = SMILES('[CH]=[C]CC(C)=C[O]'),
    E0 = (467.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.728348,'amu*angstrom^2'), symmetry=1, barrier=(16.7461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.728684,'amu*angstrom^2'), symmetry=1, barrier=(16.7539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.728633,'amu*angstrom^2'), symmetry=1, barrier=(16.7527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.641683,0.065208,-5.00835e-05,1.32287e-08,1.04609e-12,56332.1,27.7802], Tmin=(100,'K'), Tmax=(1005.17,'K')), NASAPolynomial(coeffs=[16.2516,0.0197723,-7.17567e-06,1.28206e-09,-8.9204e-14,52351.1,-51.8015], Tmin=(1005.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ) + radical(Cds_S)"""),
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
    label = '[CH]=CCC([CH2])=C[O](19971)',
    structure = SMILES('[CH]=CCC([CH2])=C[O]'),
    E0 = (380.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.710318,0.0576843,-1.51973e-05,-3.14893e-08,1.8959e-11,45950.7,27.3409], Tmin=(100,'K'), Tmax=(949.132,'K')), NASAPolynomial(coeffs=[19.4578,0.0146846,-4.14943e-06,7.23019e-10,-5.42697e-14,40770,-70.6703], Tmin=(949.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P) + radical(C=COJ)"""),
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
    label = '[CH]C(=C)C[C]=[CH](18821)',
    structure = SMILES('[CH]C(=C)C[C]=[CH]'),
    E0 = (905.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,351.212,351.261,351.665,351.728],'cm^-1')),
        HinderedRotor(inertia=(0.581927,'amu*angstrom^2'), symmetry=1, barrier=(51.055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.582794,'amu*angstrom^2'), symmetry=1, barrier=(51.0423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.581489,'amu*angstrom^2'), symmetry=1, barrier=(51.0473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3226.68,'J/mol'), sigma=(5.67268,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=504.00 K, Pc=40.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46847,0.0539599,-4.0095e-05,1.60924e-08,-2.72865e-12,108977,24.6512], Tmin=(100,'K'), Tmax=(1348.57,'K')), NASAPolynomial(coeffs=[10.0073,0.0286331,-1.19248e-05,2.16668e-09,-1.47115e-13,106674,-19.0976], Tmin=(1348.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(905.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC([CH2])=C[O](19973)',
    structure = SMILES('[CH]=[C]CC([CH2])=C[O]'),
    E0 = (618.802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,273.131,273.81],'cm^-1')),
        HinderedRotor(inertia=(0.425078,'amu*angstrom^2'), symmetry=1, barrier=(22.8387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.424101,'amu*angstrom^2'), symmetry=1, barrier=(22.8304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.425145,'amu*angstrom^2'), symmetry=1, barrier=(22.8309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676551,0.0622972,-4.00564e-05,-2.07874e-09,7.90601e-12,74554.1,27.9356], Tmin=(100,'K'), Tmax=(957.806,'K')), NASAPolynomial(coeffs=[18.2744,0.0141737,-4.42049e-06,7.70418e-10,-5.54612e-14,70019.4,-62.2795], Tmin=(957.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(C=COJ) + radical(Allyl_P)"""),
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
    label = '[CH]=[C]CC([CH2])=[C]O(28748)',
    structure = SMILES('[CH]=[C]CC([CH2])=[C]O'),
    E0 = (717.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,302.952],'cm^-1')),
        HinderedRotor(inertia=(0.266361,'amu*angstrom^2'), symmetry=1, barrier=(17.3873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265561,'amu*angstrom^2'), symmetry=1, barrier=(17.3848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265414,'amu*angstrom^2'), symmetry=1, barrier=(17.3864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267776,'amu*angstrom^2'), symmetry=1, barrier=(17.3882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.229545,0.0739095,-8.03635e-05,4.33507e-08,-8.98668e-12,86388.7,31.2584], Tmin=(100,'K'), Tmax=(1231.1,'K')), NASAPolynomial(coeffs=[18.4045,0.0129135,-3.6769e-06,5.4115e-10,-3.2963e-14,82060.9,-59.6061], Tmin=(1231.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(717.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]C[C]1CC1O(28749)',
    structure = SMILES('[CH]=[C]C[C]1CC1O'),
    E0 = (554.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31962,0.0465485,-2.00055e-06,-3.25897e-08,1.62362e-11,66800.2,26.8913], Tmin=(100,'K'), Tmax=(981.565,'K')), NASAPolynomial(coeffs=[14.6525,0.0212506,-7.71219e-06,1.42594e-09,-1.02976e-13,62784,-44.3093], Tmin=(981.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJ(C)CO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C[C]([CH2])C1O(28750)',
    structure = SMILES('[CH]=C1C[C]([CH2])C1O'),
    E0 = (499.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.8414,-0.0163826,0.000129105,-1.2344e-07,2.83903e-11,59709.6,-14.5198], Tmin=(100,'K'), Tmax=(1718.81,'K')), NASAPolynomial(coeffs=[74.3538,0.0233251,-6.91907e-05,1.69435e-08,-1.26034e-12,10178,-435.906], Tmin=(1718.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Isobutyl) + radical(Cds_P) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2][C]1C[C]=CC1O(28751)',
    structure = SMILES('[CH2][C]1C[C]=CC1O'),
    E0 = (431.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85267,0.0300579,4.60018e-05,-8.29795e-08,3.44891e-11,51994.2,25.3611], Tmin=(100,'K'), Tmax=(951.642,'K')), NASAPolynomial(coeffs=[14.4322,0.0201978,-6.2577e-06,1.12822e-09,-8.412e-14,47652.2,-44.9381], Tmin=(951.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]=C=CC(C)=CO(28752)',
    structure = SMILES('[CH]=C=CC(C)=CO'),
    E0 = (142.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11335,0.0842604,-8.80621e-05,4.3593e-08,-7.88307e-12,17358.4,28.3095], Tmin=(100,'K'), Tmax=(1626.52,'K')), NASAPolynomial(coeffs=[22.5557,0.00647779,1.72246e-06,-6.06411e-10,4.77607e-14,12248.1,-89.4337], Tmin=(1626.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]C[C]([CH2])C[O](19981)',
    structure = SMILES('[CH]=[C]C[C]([CH2])C[O]'),
    E0 = (835.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3120,650,792.5,1650,360,370,350,180,935.855,939.193],'cm^-1')),
        HinderedRotor(inertia=(0.00522917,'amu*angstrom^2'), symmetry=1, barrier=(3.27661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143586,'amu*angstrom^2'), symmetry=1, barrier=(3.30133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142773,'amu*angstrom^2'), symmetry=1, barrier=(3.28263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142471,'amu*angstrom^2'), symmetry=1, barrier=(3.27569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18674,0.0675353,-8.74487e-05,7.39128e-08,-2.5887e-11,100565,30.8037], Tmin=(100,'K'), Tmax=(823.865,'K')), NASAPolynomial(coeffs=[5.24809,0.0383929,-1.72317e-05,3.2096e-09,-2.19196e-13,100216,13.9379], Tmin=(823.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(835.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P) + radical(CCOJ) + radical(CCJ(C)CO) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CC([CH2])[CH]O(19983)',
    structure = SMILES('[CH][C]=CC([CH2])[CH]O'),
    E0 = (749.779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44633,0.081484,-0.000101725,7.50314e-08,-2.26152e-11,90302.4,32.2579], Tmin=(100,'K'), Tmax=(860.666,'K')), NASAPolynomial(coeffs=[9.71708,0.0348581,-1.42945e-05,2.53021e-09,-1.67573e-13,88837.7,-10.3155], Tmin=(860.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(749.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH]C(=C)C([CH2])=CO(20097)',
    structure = SMILES('[CH]C(=C)C([CH2])=CO'),
    E0 = (332.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.99495,'amu*angstrom^2'), symmetry=1, barrier=(45.8679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99695,'amu*angstrom^2'), symmetry=1, barrier=(45.9139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.998,'amu*angstrom^2'), symmetry=1, barrier=(45.9379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99719,'amu*angstrom^2'), symmetry=1, barrier=(45.9193,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00410301,0.0696419,-1.74234e-05,-4.63828e-08,2.87225e-11,40107.5,24.8186], Tmin=(100,'K'), Tmax=(913.308,'K')), NASAPolynomial(coeffs=[24.6662,0.0104835,-5.58477e-07,-8.1577e-11,4.68258e-15,33562.1,-103.127], Tmin=(913.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC[C]=CO(20027)',
    structure = SMILES('[CH]=[C]CC[C]=CO'),
    E0 = (577.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.743931,'amu*angstrom^2'), symmetry=1, barrier=(17.1044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745038,'amu*angstrom^2'), symmetry=1, barrier=(17.1299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.742513,'amu*angstrom^2'), symmetry=1, barrier=(17.0718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743422,'amu*angstrom^2'), symmetry=1, barrier=(17.0927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3854.44,'J/mol'), sigma=(6.42147,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.05 K, Pc=33.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0439716,0.0768172,-7.99886e-05,4.1189e-08,-8.11398e-12,69635.6,30.7676], Tmin=(100,'K'), Tmax=(1337.11,'K')), NASAPolynomial(coeffs=[19.7639,0.0130433,-3.37705e-06,4.64312e-10,-2.71878e-14,64742.5,-69.0382], Tmin=(1337.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C1CC(=CO)C1(28753)',
    structure = SMILES('[CH]=C1CC(=CO)C1'),
    E0 = (234.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21889,0.0337151,6.91893e-05,-1.31622e-07,5.79839e-11,28322.2,23.2683], Tmin=(100,'K'), Tmax=(923.677,'K')), NASAPolynomial(coeffs=[24.5285,0.0036768,2.82515e-06,-6.17717e-10,3.32135e-14,20991.3,-103.71], Tmin=(923.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC([CH2])C=O(18136)',
    structure = SMILES('[CH]=[C]CC([CH2])C=O'),
    E0 = (542.806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,3000,3100,440,815,1455,1000,301.069],'cm^-1')),
        HinderedRotor(inertia=(0.111773,'amu*angstrom^2'), symmetry=1, barrier=(7.12128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112059,'amu*angstrom^2'), symmetry=1, barrier=(7.12534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110752,'amu*angstrom^2'), symmetry=1, barrier=(7.14336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205849,'amu*angstrom^2'), symmetry=1, barrier=(13.2408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3598.66,'J/mol'), sigma=(6.12618,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.10 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.747022,0.077154,-0.00010676,8.74946e-08,-2.93851e-11,65396.3,29.3671], Tmin=(100,'K'), Tmax=(794.046,'K')), NASAPolynomial(coeffs=[8.24586,0.0342149,-1.58907e-05,3.01277e-09,-2.0804e-13,64368.2,-4.05612], Tmin=(794.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([CH2])=CO(17684)',
    structure = SMILES('[CH2]C([CH2])=CO'),
    E0 = (61.2744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.31967,'amu*angstrom^2'), symmetry=1, barrier=(30.3418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32161,'amu*angstrom^2'), symmetry=1, barrier=(30.3865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32215,'amu*angstrom^2'), symmetry=1, barrier=(30.3989,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44673,0.0388097,1.68219e-05,-6.63718e-08,3.33598e-11,7477.77,17.5086], Tmin=(100,'K'), Tmax=(909.843,'K')), NASAPolynomial(coeffs=[20.4722,0.00116961,3.03528e-06,-6.98811e-10,4.52843e-14,2111.65,-82.9443], Tmin=(909.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.2744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[C]=[CH](18830)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([562.459,1392.74,3112.83],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86574,0.0023643,1.31489e-06,-2.64796e-09,1.00932e-12,101918,6.62295], Tmin=(100,'K'), Tmax=(1016.96,'K')), NASAPolynomial(coeffs=[4.17915,0.00251347,-9.43473e-07,1.6872e-10,-1.15831e-14,101782,4.75427], Tmin=(1016.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
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
    label = '[CH]=[C]C[C]=CO(28205)',
    structure = SMILES('[CH]=[C]C[C]=CO'),
    E0 = (602.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.862394,'amu*angstrom^2'), symmetry=1, barrier=(19.8281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86231,'amu*angstrom^2'), symmetry=1, barrier=(19.8262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863234,'amu*angstrom^2'), symmetry=1, barrier=(19.8475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.552409,0.0625021,-6.76301e-05,3.49309e-08,-6.72942e-12,72628,26.4093], Tmin=(100,'K'), Tmax=(1449.98,'K')), NASAPolynomial(coeffs=[17.9977,0.0058566,-2.16834e-07,-1.16428e-10,1.19744e-14,68464.6,-61.1478], Tmin=(1449.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=CO)C[C]=[CH](28754)',
    structure = SMILES('[CH]C(=CO)C[C]=[CH]'),
    E0 = (696.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.04699,'amu*angstrom^2'), symmetry=1, barrier=(47.0643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04695,'amu*angstrom^2'), symmetry=1, barrier=(47.0634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04673,'amu*angstrom^2'), symmetry=1, barrier=(47.0583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04711,'amu*angstrom^2'), symmetry=1, barrier=(47.067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.399971,0.0793833,-7.78988e-05,3.77149e-08,-6.94236e-12,83945.8,31.2269], Tmin=(100,'K'), Tmax=(1488.65,'K')), NASAPolynomial(coeffs=[20.7528,0.0139206,-3.24576e-06,3.90526e-10,-2.05507e-14,78603.7,-76.0288], Tmin=(1488.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=[C]CC([CH2])=CO(28755)',
    structure = SMILES('[C]=[C]CC([CH2])=CO'),
    E0 = (788.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.888868,'amu*angstrom^2'), symmetry=1, barrier=(20.4368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.88878,'amu*angstrom^2'), symmetry=1, barrier=(20.4348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889128,'amu*angstrom^2'), symmetry=1, barrier=(20.4428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888952,'amu*angstrom^2'), symmetry=1, barrier=(20.4388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.198936,0.0779757,-8.42149e-05,4.36259e-08,-8.5128e-12,94979.5,29.6455], Tmin=(100,'K'), Tmax=(1393.22,'K')), NASAPolynomial(coeffs=[21.3137,0.00863934,-1.41145e-06,1.02434e-10,-2.88199e-15,89720.1,-78.6369], Tmin=(1393.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(788.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC(=C)C=O(19965)',
    structure = SMILES('[CH]=[C]CC(=C)C=O'),
    E0 = (442.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1685,370,310.195],'cm^-1')),
        HinderedRotor(inertia=(0.182925,'amu*angstrom^2'), symmetry=1, barrier=(12.488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182881,'amu*angstrom^2'), symmetry=1, barrier=(12.488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182856,'amu*angstrom^2'), symmetry=1, barrier=(12.4879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68419,0.0569697,-5.07712e-05,2.53407e-08,-5.58423e-12,53344.7,24.6538], Tmin=(100,'K'), Tmax=(1027.47,'K')), NASAPolynomial(coeffs=[7.63504,0.0338027,-1.69497e-05,3.39579e-09,-2.44662e-13,52121.9,-4.2169], Tmin=(1027.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC(=C)C[O](20704)',
    structure = SMILES('[CH]=[C]CC(=C)C[O]'),
    E0 = (607.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,322.176,322.194,1847.02],'cm^-1')),
        HinderedRotor(inertia=(0.0921395,'amu*angstrom^2'), symmetry=1, barrier=(6.78756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092148,'amu*angstrom^2'), symmetry=1, barrier=(6.78771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0921641,'amu*angstrom^2'), symmetry=1, barrier=(6.78791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25205,0.0662074,-8.51446e-05,7.31475e-08,-2.65782e-11,73122.9,28.9367], Tmin=(100,'K'), Tmax=(782.093,'K')), NASAPolynomial(coeffs=[4.99908,0.0392373,-1.84467e-05,3.5316e-09,-2.45749e-13,72775.6,13.3066], Tmin=(782.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C=C([CH2])CO(28756)',
    structure = SMILES('[CH]=[C]C=C([CH2])CO'),
    E0 = (429.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716025,0.068894,-6.93782e-05,3.84628e-08,-8.56651e-12,51820.4,27.6562], Tmin=(100,'K'), Tmax=(1091.79,'K')), NASAPolynomial(coeffs=[12.7521,0.0247965,-8.79214e-06,1.46723e-09,-9.50619e-14,49192.2,-31.4678], Tmin=(1091.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C]CC(=[CH])CO(28757)',
    structure = SMILES('[CH]=[C]CC(=[CH])CO'),
    E0 = (628.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3115,3125,620,680,785,800,1600,1700,350,440,435,1725,342.174],'cm^-1')),
        HinderedRotor(inertia=(0.102426,'amu*angstrom^2'), symmetry=1, barrier=(8.53554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102828,'amu*angstrom^2'), symmetry=1, barrier=(8.54352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102319,'amu*angstrom^2'), symmetry=1, barrier=(8.53981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10253,'amu*angstrom^2'), symmetry=1, barrier=(8.5465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.935567,0.0709073,-8.31421e-05,5.56165e-08,-1.53207e-11,75709.5,29.6251], Tmin=(100,'K'), Tmax=(876.74,'K')), NASAPolynomial(coeffs=[9.88472,0.0300769,-1.32841e-05,2.49554e-09,-1.72964e-13,74140.4,-12.3719], Tmin=(876.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(628.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=CO)CC=[CH](28758)',
    structure = SMILES('[CH]C(=CO)CC=[CH]'),
    E0 = (458.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.304711,0.0670036,-2.65109e-05,-2.52777e-08,1.81152e-11,55313,28.2151], Tmin=(100,'K'), Tmax=(931.669,'K')), NASAPolynomial(coeffs=[20.3413,0.0171757,-4.56445e-06,7.19199e-10,-5.05546e-14,50008.5,-75.4627], Tmin=(931.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    label = '[CH]=C1C[C]([CH]O)C1(28759)',
    structure = SMILES('[CH]=C1C[C]([CH]O)C1'),
    E0 = (487.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.4675,-0.014745,0.000122797,-1.17784e-07,2.67905e-11,58253.7,-19.6005], Tmin=(100,'K'), Tmax=(1751.39,'K')), NASAPolynomial(coeffs=[80.2115,0.0180778,-6.78935e-05,1.66882e-08,-1.23835e-12,5060.57,-472.891], Tmin=(1751.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CCJ(C)CO) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = 'O[CH][C]1C[C]=CC1(28760)',
    structure = SMILES('O[CH][C]1C[C]=CC1'),
    E0 = (419.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54621,0.0421626,5.2753e-06,-3.53286e-08,1.59516e-11,50580.2,24.3588], Tmin=(100,'K'), Tmax=(1000.86,'K')), NASAPolynomial(coeffs=[12.7509,0.0246973,-9.48677e-06,1.77271e-09,-1.27054e-13,46969.2,-36.5417], Tmin=(1000.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CCsJOH) + radical(cyclopentene-vinyl) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C#CC=C([CH2])CO(28761)',
    structure = SMILES('C#CC=C([CH2])CO'),
    E0 = (163.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995649,0.0639752,-5.9105e-05,2.99442e-08,-6.20383e-12,19823.4,23.6409], Tmin=(100,'K'), Tmax=(1154.69,'K')), NASAPolynomial(coeffs=[11.797,0.026557,-1.04957e-05,1.8787e-09,-1.2727e-13,17329.1,-30.0226], Tmin=(1154.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]=CCC(=C)C=O(19979)',
    structure = SMILES('[CH]=CCC(=C)C=O'),
    E0 = (205.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43714,0.0558737,-3.88167e-05,1.29687e-08,-1.76035e-12,24752.5,25.0478], Tmin=(100,'K'), Tmax=(1659.23,'K')), NASAPolynomial(coeffs=[13.5328,0.0267139,-1.24552e-05,2.37677e-09,-1.64441e-13,20738.6,-39.4316], Tmin=(1659.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC(O)[C]=C(20061)',
    structure = SMILES('[CH]=[C]CC(O)[C]=C'),
    E0 = (613.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1670,1700,300,440,287.907,288.438],'cm^-1')),
        HinderedRotor(inertia=(0.00210941,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198883,'amu*angstrom^2'), symmetry=1, barrier=(11.8025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203226,'amu*angstrom^2'), symmetry=1, barrier=(11.7828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197077,'amu*angstrom^2'), symmetry=1, barrier=(11.7571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.5,'J/mol'), sigma=(6.3848,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.16 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01939,0.0670476,-6.81196e-05,3.76161e-08,-8.52212e-12,73913.9,30.3635], Tmin=(100,'K'), Tmax=(1055.49,'K')), NASAPolynomial(coeffs=[11.3942,0.0277299,-1.22434e-05,2.32358e-09,-1.62811e-13,71723.8,-20.2492], Tmin=(1055.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC(=C)C1O(28762)',
    structure = SMILES('[CH]=C1CC(=C)C1O'),
    E0 = (265.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56704,0.0348318,3.96223e-05,-8.01087e-08,3.389e-11,32024.6,24.5431], Tmin=(100,'K'), Tmax=(964.006,'K')), NASAPolynomial(coeffs=[16.8172,0.0178125,-5.87511e-06,1.13337e-09,-8.79882e-14,26934.9,-59.6197], Tmin=(964.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P)"""),
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
    label = '[CH]=[C]C[C]=C(17763)',
    structure = SMILES('[CH]=[C]C[C]=C'),
    E0 = (811.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,250.834],'cm^-1')),
        HinderedRotor(inertia=(0.247665,'amu*angstrom^2'), symmetry=1, barrier=(11.0543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247731,'amu*angstrom^2'), symmetry=1, barrier=(11.0546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3393,0.0378852,-3.20383e-05,1.55074e-08,-3.23394e-12,97663.1,20.1382], Tmin=(100,'K'), Tmax=(1106.63,'K')), NASAPolynomial(coeffs=[7.00928,0.0210051,-9.15799e-06,1.72363e-09,-1.2002e-13,96629.5,-2.86492], Tmin=(1106.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(811.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(Cds_S)"""),
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
    E0 = (477.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (595.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (557.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (504.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (542.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (620.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (594.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (679.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (732.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (521.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (871.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (636.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (582.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (541.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (628.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (933.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (830.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (801.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (758.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (929.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (708.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (603.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (509.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (540.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (858.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (772.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (722.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (747.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (485.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (621.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (942.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1018.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (908.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1000.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (675.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (719.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (653.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (773.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (871.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (668.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (603.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (509.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (540.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (502.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (708.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (485.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (1051.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['C=C=CO(12571)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=C1CC1([CH2])[CH]O(28742)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(117.712,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 117.5 to 117.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH2]C1([CH]O)C=[C]C1(28743)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(80.6246,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R5;doublebond_intra_HNd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 79.0 to 80.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC=C([CH2])[CH]O(27790)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(197.535,'m^3/(mol*s)'), n=1.6075, Ea=(11.7084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][C]=CO(18753)', 'C#C[CH2](17441)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00755369,'m^3/(mol*s)'), n=2.43138, Ea=(27.0165,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C=CO(12571)', '[CH][C]=C(18825)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]C=CC([CH2])=CO(28744)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18474e+10,'s^-1'), n=0.962222, Ea=(117.477,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['C#C[CH][C](C)[CH]O(28745)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=[C]CC(C)=[C]O(28746)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=CCC([CH2])=[C]O(28747)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out;Cd_H_out_singleNd]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=[C]CC(C)=C[O](19967)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(493835,'s^-1'), n=1.76395, Ea=(159.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_2H;XH_out] for rate rule [R4H_SDS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH2]C(=[C]O)C[C]=C(28593)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=CCC([CH2])=C[O](19971)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(313415,'s^-1'), n=1.7968, Ea=(63.8264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_RSMS;Cd_rad_out;XH_out] + [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH2]C(=C[O])C[C]=C(18133)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OH(D)(132)', '[CH]C(=C)C[C]=[CH](18821)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_pri_rad] for rate rule [Cd_rad;O_pri_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH]=[C]CC([CH2])=C[O](19973)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [O_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=CO(18753)', '[CH][C]=C(18825)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH]=[C]C=C([CH2])[CH]O(21903)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[CH]=[C]CC([CH2])=[C]O(28748)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=[C]C[C]1CC1O(28749)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=C1C[C]([CH2])C1O(28750)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH2][C]1C[C]=CC1O(28751)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=C=CC(C)=CO(28752)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]C[C]([CH2])C[O](19981)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH][C]=CC([CH2])[CH]O(19983)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]C(=C)C([CH2])=CO(20097)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]CC[C]=CO(20027)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=C1CC(=CO)C1(28753)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=[C]CC([CH2])C=O(18136)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([CH2])=CO(17684)', '[C]=[CH](18830)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(T)(28)', '[CH]=[C]C[C]=CO(28205)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(8)', '[CH]C(=CO)C[C]=[CH](28754)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(8)', '[C]=[C]CC([CH2])=CO(28755)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(8)', '[CH]=[C]CC(=C)C=O(19965)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2834 used for Od_CO-CdH;HJ
Exact match found for rate rule [Od_CO-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=[C]CC(=C)C[O](20704)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=[C]C=C([CH2])CO(28756)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.00703183,'s^-1'), n=4.63833, Ea=(175.937,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeO;XH_out] + [R3H_SS_2Cd;C_rad_out_1H;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C]CC(=[CH])CO(28757)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]C(=CO)CC=[CH](28758)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.80239e+12,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]C(=CO)C[C]=C(28598)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=C1C[C]([CH]O)C1(28759)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['O[CH][C]1C[C]=CC1(28760)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['C#CC=C([CH2])CO(28761)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=CCC(=C)C=O(19979)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]CC(O)[C]=C(20061)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=[C]CC([CH2])=CO(19986)'],
    products = ['[CH]=C1CC(=C)C1O(28762)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]O(5471)', '[CH]=[C]C[C]=C(17763)'],
    products = ['[CH]=[C]CC([CH2])=CO(19986)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

network(
    label = '5238',
    isomers = [
        '[CH]=[C]CC([CH2])=CO(19986)',
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
    label = '5238',
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

