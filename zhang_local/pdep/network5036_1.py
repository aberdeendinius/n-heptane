species(
    label = '[CH2]C([O])C[C]=CO(14779)',
    structure = SMILES('[CH2]C([O])C[C]=CO'),
    E0 = (258.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,285.171,285.312,285.729],'cm^-1')),
        HinderedRotor(inertia=(0.291222,'amu*angstrom^2'), symmetry=1, barrier=(16.8239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291172,'amu*angstrom^2'), symmetry=1, barrier=(16.8248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291385,'amu*angstrom^2'), symmetry=1, barrier=(16.8249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291982,'amu*angstrom^2'), symmetry=1, barrier=(16.8191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.435558,0.0818472,-8.73307e-05,4.49361e-08,-8.705e-12,31230.6,31.3784], Tmin=(100,'K'), Tmax=(1412.54,'K')), NASAPolynomial(coeffs=[21.9992,0.00966968,-1.50175e-06,9.41924e-11,-1.57953e-15,25755.3,-81.5518], Tmin=(1412.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = 'C=C[O](594)',
    structure = SMILES('C=C[O]'),
    E0 = (-25.1807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34719,0.00128739,5.39982e-05,-7.84138e-08,3.24083e-11,-2992.85,8.97297], Tmin=(100,'K'), Tmax=(914.213,'K')), NASAPolynomial(coeffs=[11.726,-0.0014735,2.90737e-06,-5.96989e-10,3.70275e-14,-5941.49,-38.4465], Tmin=(914.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.1807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
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
    label = 'C=C([O])C[C]=CO(14467)',
    structure = SMILES('C=C([O])C[C]=CO'),
    E0 = (41.0451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.749384,'amu*angstrom^2'), symmetry=1, barrier=(17.2298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750192,'amu*angstrom^2'), symmetry=1, barrier=(17.2484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75198,'amu*angstrom^2'), symmetry=1, barrier=(17.2895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4284.46,'J/mol'), sigma=(6.81655,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.22 K, Pc=30.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.56428,0.0758759,-7.90326e-05,3.87517e-08,-6.96524e-12,5122.04,31.1132], Tmin=(100,'K'), Tmax=(1620.33,'K')), NASAPolynomial(coeffs=[21.4332,0.00536805,1.23983e-06,-4.47392e-10,3.50195e-14,120.562,-79.0647], Tmin=(1620.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.0451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C=C=CO(23314)',
    structure = SMILES('[CH2]C([O])C=C=CO'),
    E0 = (184.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0094,'amu*angstrom^2'), symmetry=1, barrier=(23.208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00858,'amu*angstrom^2'), symmetry=1, barrier=(23.1893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01006,'amu*angstrom^2'), symmetry=1, barrier=(23.2232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.783812,0.0797207,-8.30818e-05,4.01636e-08,-7.13092e-12,22373.3,31.6978], Tmin=(100,'K'), Tmax=(1624.42,'K')), NASAPolynomial(coeffs=[23.4588,0.00425884,1.15907e-06,-3.84135e-10,2.89781e-14,16577.4,-90.6173], Tmin=(1624.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])CC#CO(28661)',
    structure = SMILES('[CH2]C([O])CC#CO'),
    E0 = (253.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2100,2250,500,550,1380,1390,370,380,2900,435,301.977,301.989,301.995],'cm^-1')),
        HinderedRotor(inertia=(0.00184838,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00184872,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.334018,'amu*angstrom^2'), symmetry=1, barrier=(21.6171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1572,'amu*angstrom^2'), symmetry=1, barrier=(74.8835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.971913,0.0655508,-6.85246e-05,3.84511e-08,-8.67575e-12,30636.9,27.1735], Tmin=(100,'K'), Tmax=(1073.82,'K')), NASAPolynomial(coeffs=[12.3554,0.023147,-9.29147e-06,1.67703e-09,-1.14239e-13,28192.1,-28.5562], Tmin=(1073.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CtH) + group(Cs-CsCsOsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = 'O=CC[C]=CO(27745)',
    structure = SMILES('O=CC[C]=CO'),
    E0 = (-64.0896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.935575,'amu*angstrom^2'), symmetry=1, barrier=(21.5107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.93612,'amu*angstrom^2'), symmetry=1, barrier=(21.5232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.935243,'amu*angstrom^2'), symmetry=1, barrier=(21.5031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49226,0.0417242,-4.4337e-06,-3.09691e-08,1.61113e-11,-7605.93,20.9994], Tmin=(100,'K'), Tmax=(988.339,'K')), NASAPolynomial(coeffs=[17.1066,0.00994677,-3.8861e-06,8.23847e-10,-6.62225e-14,-12226.8,-61.91], Tmin=(988.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.0896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(Cds_S)"""),
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
    label = 'C=CC[C]=CO(27709)',
    structure = SMILES('C=CC[C]=CO'),
    E0 = (117.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.86248,'amu*angstrom^2'), symmetry=1, barrier=(19.8301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86521,'amu*angstrom^2'), symmetry=1, barrier=(19.8929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.865876,'amu*angstrom^2'), symmetry=1, barrier=(19.9082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17072,0.0483781,-5.30116e-06,-3.87313e-08,2.16764e-11,14282.6,23.3117], Tmin=(100,'K'), Tmax=(927.027,'K')), NASAPolynomial(coeffs=[18.1948,0.0105692,-1.80457e-06,2.3518e-10,-1.86167e-14,9594.49,-65.7914], Tmin=(927.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = '[CH2][CH][O](719)',
    structure = SMILES('[CH2][CH][O]'),
    E0 = (361.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1878.99],'cm^-1')),
        HinderedRotor(inertia=(0.232981,'amu*angstrom^2'), symmetry=1, barrier=(5.35669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03639,0.0272039,-5.17476e-05,5.40082e-08,-2.05139e-11,43449.8,12.3205], Tmin=(100,'K'), Tmax=(879.689,'K')), NASAPolynomial(coeffs=[2.12305,0.0164211,-7.89343e-06,1.47303e-09,-9.88046e-14,44188.4,19.8945], Tmin=(879.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CJCO) + radical(CCOJ)"""),
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
    label = 'C#CCC([CH2])[O](17828)',
    structure = SMILES('C#CCC([CH2])[O]'),
    E0 = (395.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3000,3100,440,815,1455,1000,371.211,372.118],'cm^-1')),
        HinderedRotor(inertia=(0.0681491,'amu*angstrom^2'), symmetry=1, barrier=(6.66411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21733,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21846,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1544,0.0586267,-5.99396e-05,3.28582e-08,-7.14731e-12,47647.1,23.1773], Tmin=(100,'K'), Tmax=(1123.69,'K')), NASAPolynomial(coeffs=[12.4442,0.0184382,-6.29232e-06,1.03001e-09,-6.60987e-14,45109.9,-32.6059], Tmin=(1123.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](O)C[C]=CO(28362)',
    structure = SMILES('[CH2][C](O)C[C]=CO'),
    E0 = (204.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,270.153,270.156],'cm^-1')),
        HinderedRotor(inertia=(0.266841,'amu*angstrom^2'), symmetry=1, barrier=(13.8195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266853,'amu*angstrom^2'), symmetry=1, barrier=(13.8194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26683,'amu*angstrom^2'), symmetry=1, barrier=(13.8195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26683,'amu*angstrom^2'), symmetry=1, barrier=(13.8196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26682,'amu*angstrom^2'), symmetry=1, barrier=(13.8195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.156889,0.0875983,-0.000106099,6.11056e-08,-1.27597e-11,24748.1,29.9923], Tmin=(100,'K'), Tmax=(911.222,'K')), NASAPolynomial(coeffs=[19.2903,0.0139482,-4.14865e-06,6.28611e-10,-3.91032e-14,20717.6,-64.6902], Tmin=(911.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CJCO) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = 'C[C]([O])C[C]=CO(28164)',
    structure = SMILES('C[C]([O])C[C]=CO'),
    E0 = (223.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,231.865,232.032,232.077],'cm^-1')),
        HinderedRotor(inertia=(0.392217,'amu*angstrom^2'), symmetry=1, barrier=(15.0267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.391866,'amu*angstrom^2'), symmetry=1, barrier=(15.0332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393631,'amu*angstrom^2'), symmetry=1, barrier=(15.0368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.392448,'amu*angstrom^2'), symmetry=1, barrier=(15.0322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.12061,0.07888,-8.21698e-05,3.93008e-08,-6.08764e-12,26998.1,28.1691], Tmin=(100,'K'), Tmax=(953.789,'K')), NASAPolynomial(coeffs=[18.5366,0.0157197,-4.97119e-06,8.1111e-10,-5.37598e-14,22845,-63.1618], Tmin=(953.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])[CH]C=CO(14780)',
    structure = SMILES('[CH2]C([O])[CH]C=CO'),
    E0 = (137.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14256,0.0870862,-8.97305e-05,4.32837e-08,-7.71656e-12,16721.1,30.6599], Tmin=(100,'K'), Tmax=(1602.27,'K')), NASAPolynomial(coeffs=[24.9552,0.00639511,3.57512e-07,-2.52614e-10,2.07641e-14,10352.6,-101.325], Tmin=(1602.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CJCO) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C([O])CC=[C]O(14778)',
    structure = SMILES('[CH2]C([O])CC=[C]O'),
    E0 = (260.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,297.087,297.094,297.151],'cm^-1')),
        HinderedRotor(inertia=(0.00192082,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24994,'amu*angstrom^2'), symmetry=1, barrier=(15.6394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249931,'amu*angstrom^2'), symmetry=1, barrier=(15.6348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249673,'amu*angstrom^2'), symmetry=1, barrier=(15.643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0590339,0.0760536,-7.91389e-05,4.10809e-08,-8.21293e-12,31437.2,32.2568], Tmin=(100,'K'), Tmax=(1284.29,'K')), NASAPolynomial(coeffs=[19.0484,0.0143794,-4.15001e-06,6.20395e-10,-3.82049e-14,26768.4,-63.2946], Tmin=(1284.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=CJO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)[CH][C]=CO(28363)',
    structure = SMILES('[CH2]C(O)[CH][C]=CO'),
    E0 = (144.782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02282,'amu*angstrom^2'), symmetry=1, barrier=(23.5165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02304,'amu*angstrom^2'), symmetry=1, barrier=(23.5218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02298,'amu*angstrom^2'), symmetry=1, barrier=(23.5202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02289,'amu*angstrom^2'), symmetry=1, barrier=(23.5181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02301,'amu*angstrom^2'), symmetry=1, barrier=(23.521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.799666,0.089634,-0.000100589,5.35229e-08,-1.0668e-11,17599.8,30.0735], Tmin=(100,'K'), Tmax=(1373.03,'K')), NASAPolynomial(coeffs=[24.0912,0.00760905,-5.8804e-07,-7.72634e-11,1.0127e-14,11661.1,-94.6365], Tmin=(1373.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJCO) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = 'CC([O])[CH][C]=CO(28165)',
    structure = SMILES('CC([O])[CH][C]=CO'),
    E0 = (163.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18051,'amu*angstrom^2'), symmetry=1, barrier=(27.1422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17918,'amu*angstrom^2'), symmetry=1, barrier=(27.1116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17966,'amu*angstrom^2'), symmetry=1, barrier=(27.1227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18043,'amu*angstrom^2'), symmetry=1, barrier=(27.1403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.771249,0.0836128,-8.5035e-05,4.11746e-08,-7.4556e-12,19861.1,29.1609], Tmin=(100,'K'), Tmax=(1543.29,'K')), NASAPolynomial(coeffs=[23.5446,0.00893448,-1.12337e-06,3.32219e-11,1.65899e-15,13743.8,-94.2031], Tmin=(1543.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])CC=CO(10682)',
    structure = SMILES('[CH2][C]([O])CC=CO'),
    E0 = (197.013,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,274.681,274.681,274.681],'cm^-1')),
        HinderedRotor(inertia=(0.301655,'amu*angstrom^2'), symmetry=1, barrier=(16.1508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301655,'amu*angstrom^2'), symmetry=1, barrier=(16.1508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301655,'amu*angstrom^2'), symmetry=1, barrier=(16.1508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301655,'amu*angstrom^2'), symmetry=1, barrier=(16.1508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.399136,0.0841836,-9.35662e-05,5.04094e-08,-1.02831e-11,23864.3,30.1943], Tmin=(100,'K'), Tmax=(1311.05,'K')), NASAPolynomial(coeffs=[21.6802,0.0105075,-2.04962e-06,2.00955e-10,-8.69053e-15,18617.4,-80.2363], Tmin=(1311.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])CC=C[O](10963)',
    structure = SMILES('[CH2]C([O])CC=C[O]'),
    E0 = (161.847,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,419.225,419.32,419.321,419.34],'cm^-1')),
        HinderedRotor(inertia=(0.139288,'amu*angstrom^2'), symmetry=1, barrier=(17.3847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139308,'amu*angstrom^2'), symmetry=1, barrier=(17.3853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139399,'amu*angstrom^2'), symmetry=1, barrier=(17.384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4287.55,'J/mol'), sigma=(7.04051,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.70 K, Pc=27.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.566955,0.0637544,-3.65769e-05,-7.0627e-09,9.74719e-12,19600,28.7134], Tmin=(100,'K'), Tmax=(961.369,'K')), NASAPolynomial(coeffs=[18.6295,0.0161084,-5.15484e-06,9.098e-10,-6.55933e-14,14855.9,-64.3273], Tmin=(961.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C[C]=[C]O(28364)',
    structure = SMILES('[CH2]C(O)C[C]=[C]O'),
    E0 = (267.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3580,3650,1210,1345,900,1100,1670,1700,300,440,1380,1390,370,380,2900,435,286.165,286.172],'cm^-1')),
        HinderedRotor(inertia=(0.00205848,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222324,'amu*angstrom^2'), symmetry=1, barrier=(12.9201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22233,'amu*angstrom^2'), symmetry=1, barrier=(12.9201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222322,'amu*angstrom^2'), symmetry=1, barrier=(12.9201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222329,'amu*angstrom^2'), symmetry=1, barrier=(12.9201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.185973,0.0808832,-9.68895e-05,5.89536e-08,-1.39395e-11,32326,32.4649], Tmin=(100,'K'), Tmax=(1043.79,'K')), NASAPolynomial(coeffs=[16.6523,0.017781,-6.20715e-06,1.03496e-09,-6.73391e-14,28888.5,-47.6815], Tmin=(1043.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'CC([O])C[C]=[C]O(28166)',
    structure = SMILES('CC([O])C[C]=[C]O'),
    E0 = (286.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,254.167,254.219,255.701],'cm^-1')),
        HinderedRotor(inertia=(0.00258328,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311015,'amu*angstrom^2'), symmetry=1, barrier=(14.2736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31064,'amu*angstrom^2'), symmetry=1, barrier=(14.2677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308636,'amu*angstrom^2'), symmetry=1, barrier=(14.2738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.397356,0.0729116,-7.53854e-05,3.9965e-08,-8.30038e-12,34578.9,30.8808], Tmin=(100,'K'), Tmax=(1180.63,'K')), NASAPolynomial(coeffs=[16.3816,0.0187566,-6.58105e-06,1.11325e-09,-7.34615e-14,30804.6,-48.8881], Tmin=(1180.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C[C]=C[O](13903)',
    structure = SMILES('[CH2]C(O)C[C]=C[O]'),
    E0 = (169.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,301.525,301.529,301.571],'cm^-1')),
        HinderedRotor(inertia=(0.00185302,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253893,'amu*angstrom^2'), symmetry=1, barrier=(16.3886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253803,'amu*angstrom^2'), symmetry=1, barrier=(16.3887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253843,'amu*angstrom^2'), symmetry=1, barrier=(16.3885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.123842,0.0752098,-7.69979e-05,3.94101e-08,-7.80785e-12,20513.6,30.9723], Tmin=(100,'K'), Tmax=(1256.84,'K')), NASAPolynomial(coeffs=[18.7481,0.0153526,-4.86324e-06,7.77955e-10,-4.99294e-14,15878.2,-62.9529], Tmin=(1256.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CJCO) + radical(C=COJ)"""),
)

species(
    label = 'CC([O])C[C]=C[O](13735)',
    structure = SMILES('CC([O])C[C]=C[O]'),
    E0 = (188.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,390.151,390.282,390.297,390.339],'cm^-1')),
        HinderedRotor(inertia=(0.147072,'amu*angstrom^2'), symmetry=1, barrier=(15.8772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14708,'amu*angstrom^2'), symmetry=1, barrier=(15.8777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146886,'amu*angstrom^2'), symmetry=1, barrier=(15.8766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674313,0.0633024,-4.20767e-05,3.61142e-09,4.71495e-12,22751.7,28.1681], Tmin=(100,'K'), Tmax=(988.706,'K')), NASAPolynomial(coeffs=[16.7385,0.0192004,-6.859e-06,1.23346e-09,-8.69739e-14,18554.2,-54.313], Tmin=(988.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]CC([CH2])[O](17830)',
    structure = SMILES('[CH]=[C]CC([CH2])[O]'),
    E0 = (714.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3120,650,792.5,1650,3000,3100,440,815,1455,1000,358.676,359.105],'cm^-1')),
        HinderedRotor(inertia=(0.00130755,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134185,'amu*angstrom^2'), symmetry=1, barrier=(12.2634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134345,'amu*angstrom^2'), symmetry=1, barrier=(12.2604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.28011,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34582,0.0599191,-6.53843e-05,3.89316e-08,-9.4385e-12,85982.6,25.5327], Tmin=(100,'K'), Tmax=(994.707,'K')), NASAPolynomial(coeffs=[10.4569,0.0232814,-1.01362e-05,1.90406e-09,-1.325e-13,84170,-18.3748], Tmin=(994.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C[C]=C[O](14477)',
    structure = SMILES('[CH2]C([O])C[C]=C[O]'),
    E0 = (399.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,430.496,430.584,430.608,430.639],'cm^-1')),
        HinderedRotor(inertia=(0.123564,'amu*angstrom^2'), symmetry=1, barrier=(16.2584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123538,'amu*angstrom^2'), symmetry=1, barrier=(16.2577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123587,'amu*angstrom^2'), symmetry=1, barrier=(16.2589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519967,0.0685223,-6.19679e-05,2.3012e-08,-1.57304e-12,48203.9,29.3555], Tmin=(100,'K'), Tmax=(988.473,'K')), NASAPolynomial(coeffs=[17.519,0.0154745,-5.35534e-06,9.40601e-10,-6.54123e-14,44074.3,-56.3477], Tmin=(988.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]([O])C[C]=CO(28662)',
    structure = SMILES('[CH2][C]([O])C[C]=CO'),
    E0 = (434.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,284.96,285.502,286.506],'cm^-1')),
        HinderedRotor(inertia=(0.265287,'amu*angstrom^2'), symmetry=1, barrier=(15.2011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26914,'amu*angstrom^2'), symmetry=1, barrier=(15.2001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262884,'amu*angstrom^2'), symmetry=1, barrier=(15.2085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265986,'amu*angstrom^2'), symmetry=1, barrier=(15.2083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00866844,0.0837864,-0.000100879,5.70373e-08,-1.16044e-11,52449.3,29.2678], Tmin=(100,'K'), Tmax=(923.845,'K')), NASAPolynomial(coeffs=[19.2581,0.0120979,-3.52913e-06,5.33074e-10,-3.34444e-14,48388.7,-64.8666], Tmin=(923.845,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C([O])[CH][C]=CO(28318)',
    structure = SMILES('[CH2]C([O])[CH][C]=CO'),
    E0 = (375.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20264,'amu*angstrom^2'), symmetry=1, barrier=(27.651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19971,'amu*angstrom^2'), symmetry=1, barrier=(27.5838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19998,'amu*angstrom^2'), symmetry=1, barrier=(27.59,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20316,'amu*angstrom^2'), symmetry=1, barrier=(27.6629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.685856,0.0861893,-9.6478e-05,5.0651e-08,-9.91865e-12,45302.5,29.4751], Tmin=(100,'K'), Tmax=(1408.55,'K')), NASAPolynomial(coeffs=[24.0253,0.00579604,1.68506e-08,-1.70375e-10,1.56418e-14,39354.8,-94.6086], Tmin=(1408.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C([O])C[C]=[C]O(28663)',
    structure = SMILES('[CH2]C([O])C[C]=[C]O'),
    E0 = (497.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1670,1700,300,440,1380,1390,370,380,2900,435,302.854,303.699,306.816],'cm^-1')),
        HinderedRotor(inertia=(0.00178658,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221299,'amu*angstrom^2'), symmetry=1, barrier=(14.4529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222478,'amu*angstrom^2'), symmetry=1, barrier=(14.4519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219576,'amu*angstrom^2'), symmetry=1, barrier=(14.4493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335187,0.0770517,-9.15571e-05,5.46641e-08,-1.26501e-11,60027.1,31.7373], Tmin=(100,'K'), Tmax=(1067.05,'K')), NASAPolynomial(coeffs=[16.655,0.0158738,-5.55585e-06,9.32088e-10,-6.10822e-14,56544.3,-48.0554], Tmin=(1067.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(O)C[C]=CO(28359)',
    structure = SMILES('C=C(O)C[C]=CO'),
    E0 = (-96.7597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.866236,'amu*angstrom^2'), symmetry=1, barrier=(19.9165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.865733,'amu*angstrom^2'), symmetry=1, barrier=(19.9049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867001,'amu*angstrom^2'), symmetry=1, barrier=(19.9341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867041,'amu*angstrom^2'), symmetry=1, barrier=(19.935,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.177022,0.0669937,-2.73134e-05,-3.49459e-08,2.50471e-11,-11483.8,26.9374], Tmin=(100,'K'), Tmax=(906.814,'K')), NASAPolynomial(coeffs=[25.3596,0.00232446,2.88644e-06,-7.06958e-10,4.74494e-14,-17959.3,-102.612], Tmin=(906.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.7597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C[C]=CO(27723)',
    structure = SMILES('CC(=O)C[C]=CO'),
    E0 = (-118.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1685,370,3010,987.5,1337.5,450,1655,320.071,320.071],'cm^-1')),
        HinderedRotor(inertia=(0.229322,'amu*angstrom^2'), symmetry=1, barrier=(16.6711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229321,'amu*angstrom^2'), symmetry=1, barrier=(16.6711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229321,'amu*angstrom^2'), symmetry=1, barrier=(16.6711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229321,'amu*angstrom^2'), symmetry=1, barrier=(16.6711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752872,0.0579818,-2.42368e-05,-1.55056e-08,1.09893e-11,-14178.9,25.246], Tmin=(100,'K'), Tmax=(1015.67,'K')), NASAPolynomial(coeffs=[18.441,0.0173621,-7.13681e-06,1.42197e-09,-1.06624e-13,-19269.9,-67.7381], Tmin=(1015.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])CC=CO(14782)',
    structure = SMILES('C=C([O])CC=CO'),
    E0 = (-196.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.564706,0.0585676,-1.06256e-05,-4.62401e-08,2.75351e-11,-23529.6,26.5726], Tmin=(100,'K'), Tmax=(911.791,'K')), NASAPolynomial(coeffs=[22.8386,0.00613579,1.13491e-06,-3.70519e-10,2.38649e-14,-29473.8,-89.1514], Tmin=(911.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C=C=CO(28360)',
    structure = SMILES('[CH2]C(O)C=C=CO'),
    E0 = (-45.9533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.882492,'amu*angstrom^2'), symmetry=1, barrier=(20.2902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.881931,'amu*angstrom^2'), symmetry=1, barrier=(20.2773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.881258,'amu*angstrom^2'), symmetry=1, barrier=(20.2618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.88246,'amu*angstrom^2'), symmetry=1, barrier=(20.2895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.871799,0.0829031,-8.64469e-05,4.22706e-08,-7.62652e-12,-5330.67,32.2004], Tmin=(100,'K'), Tmax=(1586.76,'K')), NASAPolynomial(coeffs=[23.7542,0.005773,6.94373e-07,-3.19214e-10,2.55218e-14,-11250.9,-92.0048], Tmin=(1586.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.9533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C=C=CO(28159)',
    structure = SMILES('CC([O])C=C=CO'),
    E0 = (-27.1816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.974741,'amu*angstrom^2'), symmetry=1, barrier=(22.4112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97427,'amu*angstrom^2'), symmetry=1, barrier=(22.4004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975072,'amu*angstrom^2'), symmetry=1, barrier=(22.4188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456929,0.0620193,-2.09254e-05,-3.23381e-08,2.11657e-11,-3126.97,26.589], Tmin=(100,'K'), Tmax=(932.721,'K')), NASAPolynomial(coeffs=[22.1749,0.00935299,-1.314e-06,1.65315e-10,-1.54303e-14,-8938.79,-86.1122], Tmin=(932.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.1816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ)"""),
)

species(
    label = '[O][CH]CC[C]=CO(14848)',
    structure = SMILES('[O][CH]CC[C]=CO'),
    E0 = (234.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3025,407.5,1350,352.5,275.482,275.569,275.624,275.855],'cm^-1')),
        HinderedRotor(inertia=(0.00222267,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264218,'amu*angstrom^2'), symmetry=1, barrier=(14.242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263637,'amu*angstrom^2'), symmetry=1, barrier=(14.2439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26356,'amu*angstrom^2'), symmetry=1, barrier=(14.2455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.415023,0.0797523,-9.5976e-05,6.10278e-08,-1.5409e-11,28383,28.642], Tmin=(100,'K'), Tmax=(968.612,'K')), NASAPolynomial(coeffs=[13.8704,0.0241866,-9.92643e-06,1.80229e-09,-1.22787e-13,25776.4,-35.8432], Tmin=(968.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=CO)C([CH2])[O](14760)',
    structure = SMILES('[CH2]C(=CO)C([CH2])[O]'),
    E0 = (156.273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,281.114,281.116],'cm^-1')),
        HinderedRotor(inertia=(0.384924,'amu*angstrom^2'), symmetry=1, barrier=(21.5863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.384925,'amu*angstrom^2'), symmetry=1, barrier=(21.5863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.38493,'amu*angstrom^2'), symmetry=1, barrier=(21.5863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.384928,'amu*angstrom^2'), symmetry=1, barrier=(21.5863,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0979613,0.0691298,-3.16954e-05,-2.8261e-08,2.15979e-11,18951.3,28.6894], Tmin=(100,'K'), Tmax=(921.113,'K')), NASAPolynomial(coeffs=[24.828,0.00545309,8.10985e-07,-2.64157e-10,1.50987e-14,12540.9,-98.6534], Tmin=(921.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = 'OC=[C]CC1CO1(27720)',
    structure = SMILES('OC=[C]CC1CO1'),
    E0 = (5.84038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.98275,0.0775268,-7.71134e-05,3.80145e-08,-6.7948e-12,909.623,29.9563], Tmin=(100,'K'), Tmax=(1702.56,'K')), NASAPolynomial(coeffs=[17.2696,0.0115497,1.36151e-06,-6.81058e-10,5.7068e-14,-1958.21,-57.983], Tmin=(1702.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.84038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1CC(=CO)O1(28643)',
    structure = SMILES('[CH2]C1CC(=CO)O1'),
    E0 = (-89.9432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.343564,0.0556346,2.1944e-05,-9.98349e-08,5.27284e-11,-10662.4,20.2494], Tmin=(100,'K'), Tmax=(884.887,'K')), NASAPolynomial(coeffs=[29.3749,-0.00635916,9.66375e-06,-2.15907e-09,1.51151e-13,-18511.1,-131.577], Tmin=(884.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.9432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(CJC(C)OC)"""),
)

species(
    label = '[O]C1CC(=CO)C1(28664)',
    structure = SMILES('[O]C1CC(=CO)C1'),
    E0 = (-50.5485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.59125,-0.00372062,0.000111099,-1.14641e-07,2.7041e-11,-6385.32,-15.5996], Tmin=(100,'K'), Tmax=(1704.79,'K')), NASAPolynomial(coeffs=[76.7156,0.018211,-6.60712e-05,1.63788e-08,-1.22587e-12,-55345.6,-451.715], Tmin=(1704.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.5485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(methylenecyclobutane) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][CH]C[C]=CO(27802)',
    structure = SMILES('[CH2][CH]C[C]=CO'),
    E0 = (388.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180,763.883],'cm^-1')),
        HinderedRotor(inertia=(0.030554,'amu*angstrom^2'), symmetry=1, barrier=(12.6061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21979,'amu*angstrom^2'), symmetry=1, barrier=(5.0534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304317,'amu*angstrom^2'), symmetry=1, barrier=(12.6012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724847,'amu*angstrom^2'), symmetry=1, barrier=(16.6657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79155,0.060806,-5.77015e-05,2.83828e-08,-5.43396e-12,46842.3,28.6956], Tmin=(100,'K'), Tmax=(1364.91,'K')), NASAPolynomial(coeffs=[15.1412,0.0155473,-4.44036e-06,6.47526e-10,-3.87373e-14,43223.7,-43.9031], Tmin=(1364.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = '[O][CH]C[C]=CO(28093)',
    structure = SMILES('[O][CH]C[C]=CO'),
    E0 = (258.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3025,407.5,1350,352.5,305.208,305.211,305.239],'cm^-1')),
        HinderedRotor(inertia=(0.213679,'amu*angstrom^2'), symmetry=1, barrier=(14.1255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213686,'amu*angstrom^2'), symmetry=1, barrier=(14.1255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21367,'amu*angstrom^2'), symmetry=1, barrier=(14.1253,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03666,0.0652876,-8.33377e-05,5.46181e-08,-1.39965e-11,31221.4,24.1632], Tmin=(100,'K'), Tmax=(961.256,'K')), NASAPolynomial(coeffs=[12.9274,0.0158092,-6.1313e-06,1.07434e-09,-7.15047e-14,28935.3,-32.7332], Tmin=(961.256,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([CH2])[O](740)',
    structure = SMILES('[CH2]C([CH2])[O]'),
    E0 = (360.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.00253983,'amu*angstrom^2'), symmetry=1, barrier=(28.8374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(23.1632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06989,0.0377752,-3.46167e-05,1.65319e-08,-3.12506e-12,43436.9,17.0293], Tmin=(100,'K'), Tmax=(1286.95,'K')), NASAPolynomial(coeffs=[10.5217,0.0115058,-3.99816e-06,6.70693e-10,-4.38762e-14,41261.6,-25.8777], Tmin=(1286.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJCO) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[C]=CO(27807)',
    structure = SMILES('[C]=CO'),
    E0 = (391.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.23325,'amu*angstrom^2'), symmetry=1, barrier=(28.3547,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88904,0.0147275,1.32236e-05,-4.12494e-08,2.11475e-11,47130.8,9.58665], Tmin=(100,'K'), Tmax=(884.362,'K')), NASAPolynomial(coeffs=[13.839,-0.00784764,5.80008e-06,-1.19218e-09,8.19757e-14,44140.2,-47.8537], Tmin=(884.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]C([O])C[C]=CO(28665)',
    structure = SMILES('[CH]C([O])C[C]=CO'),
    E0 = (494.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.32504,0.0795946,-8.67135e-05,4.50039e-08,-8.75355e-12,59686,31.3516], Tmin=(100,'K'), Tmax=(1413.9,'K')), NASAPolynomial(coeffs=[21.932,0.00761014,-7.7812e-07,-2.64676e-11,6.1563e-15,54293.6,-80.5467], Tmin=(1413.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCJ2_triplet) + radical(Cds_S) + radical(CC(C)OJ)"""),
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
    E0 = (258.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (274.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (404.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (480.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (373.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (360.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (258.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (351.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (423.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (372.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (399.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (415.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (492.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (333.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (376.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (400.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (408.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (382.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (319.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (331.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (331.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (742.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (611.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (547.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (646.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (586.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (709.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (281.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (281.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (347.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (336.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (336.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (415.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (428.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (263.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (266.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (266.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (795.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (674.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (786.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (706.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['C=C[O](594)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', 'C=C([O])C[C]=CO(14467)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C([O])C=C=CO(23314)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2]C([O])CC#CO(28661)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(28)', 'O=CC[C]=CO(27745)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(T)(63)', 'C=CC[C]=CO(27709)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2803 used for Cds-CsH_Cds-HH;O_atom_triplet
Exact match found for rate rule [Cds-CsH_Cds-HH;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C[O](594)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(96.7353,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 93.7 to 96.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH][O](719)', 'C=C=CO(12571)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(D)(132)', 'C#CCC([CH2])[O](17828)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.508e+07,'cm^3/(mol*s)'), n=1.628, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 211 used for Ct-H_Ct-Cs;OJ_pri
Exact match found for rate rule [Ct-H_Ct-Cs;OJ_pri]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2][C](O)C[C]=CO(28362)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['C[C]([O])C[C]=CO(28164)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C([O])[CH]C=CO(14780)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C([O])CC=[C]O(14778)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C(O)[CH][C]=CO(28363)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['CC([O])[CH][C]=CO(28165)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2][C]([O])CC=CO(10682)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C([O])CC=C[O](10963)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C(O)C[C]=[C]O(28364)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(49.9162,'s^-1'), n=2.87175, Ea=(124.402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cd_H_out_single] for rate rule [R5HJ_3;O_rad_out;Cd_H_out_singleNd]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC([O])C[C]=[C]O(28166)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C(O)C[C]=C[O](13903)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.74437e+06,'s^-1'), n=0.972854, Ea=(72.9565,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;O_rad_out;XH_out] for rate rule [R6HJ_3;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['CC([O])C[C]=C[O](13735)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['OH(D)(132)', '[CH]=[C]CC([CH2])[O](17830)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2]C([O])C[C]=C[O](14477)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH][O](719)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH2][C]([O])C[C]=CO(28662)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH2]C([O])[CH][C]=CO(28318)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', '[CH2]C([O])C[C]=[C]O(28663)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['C=C(O)C[C]=CO(28359)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['CC(=O)C[C]=CO(27723)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['C=C([O])CC=CO(14782)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C(O)C=C=CO(28360)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['CC([O])C=C=CO(28159)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[O][CH]CC[C]=CO(14848)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C(=CO)C([CH2])[O](14760)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['OC=[C]CC1CO1(27720)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[CH2]C1CC(=CO)O1(28643)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([O])C[C]=CO(14779)'],
    products = ['[O]C1CC(=CO)C1(28664)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O(T)(63)', '[CH2][CH]C[C]=CO(27802)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(T)(28)', '[O][CH]C[C]=CO(28093)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([CH2])[O](740)', '[C]=CO(27807)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['H(8)', '[CH]C([O])C[C]=CO(28665)'],
    products = ['[CH2]C([O])C[C]=CO(14779)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '5036',
    isomers = [
        '[CH2]C([O])C[C]=CO(14779)',
    ],
    reactants = [
        ('C=C[O](594)', 'C=C=CO(12571)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5036',
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

