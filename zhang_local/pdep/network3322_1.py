species(
    label = 'C[CH]CCC=C[O](12588)',
    structure = SMILES('C[CH]CCC=C[O]'),
    E0 = (63.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,365.199,1122.26,1122.32],'cm^-1')),
        HinderedRotor(inertia=(0.740781,'amu*angstrom^2'), symmetry=1, barrier=(17.032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740792,'amu*angstrom^2'), symmetry=1, barrier=(17.0323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00290702,'amu*angstrom^2'), symmetry=1, barrier=(17.0318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740774,'amu*angstrom^2'), symmetry=1, barrier=(17.0319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769553,0.0607664,-2.42806e-05,-9.80463e-09,7.66759e-12,7706.48,30.1673], Tmin=(100,'K'), Tmax=(1023.1,'K')), NASAPolynomial(coeffs=[13.6686,0.0313133,-1.18552e-05,2.14014e-09,-1.48348e-13,3969.12,-37.7236], Tmin=(1023.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJC) + radical(C=COJ)"""),
)

species(
    label = 'C=CC(42)',
    structure = SMILES('C=CC'),
    E0 = (6.12372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
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
    label = 'CC1CCC1[CH][O](13149)',
    structure = SMILES('CC1CCC1[CH][O]'),
    E0 = (194.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24715,0.0493173,1.71467e-06,-2.89211e-08,1.21705e-11,23490.1,25.4072], Tmin=(100,'K'), Tmax=(1071.52,'K')), NASAPolynomial(coeffs=[11.1004,0.0369485,-1.51465e-05,2.83251e-09,-1.99001e-13,19977,-29.3498], Tmin=(1071.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = 'CC=CCC=C[O](13150)',
    structure = SMILES('CC=CCC=C[O]'),
    E0 = (-14.6059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,206.633,206.643,206.651],'cm^-1')),
        HinderedRotor(inertia=(0.550711,'amu*angstrom^2'), symmetry=1, barrier=(16.6874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550706,'amu*angstrom^2'), symmetry=1, barrier=(16.6875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550732,'amu*angstrom^2'), symmetry=1, barrier=(16.6875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.935374,0.0527799,9.18183e-07,-4.14809e-08,2.03211e-11,-1633.24,26.8986], Tmin=(100,'K'), Tmax=(974.129,'K')), NASAPolynomial(coeffs=[16.347,0.0242059,-8.52961e-06,1.56265e-09,-1.12858e-13,-6282.67,-55.5025], Tmin=(974.129,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.6059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'C=CCCC=C[O](12800)',
    structure = SMILES('C=CCCC=C[O]'),
    E0 = (-3.63004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,330.969,332.178,332.999,333.429],'cm^-1')),
        HinderedRotor(inertia=(0.22817,'amu*angstrom^2'), symmetry=1, barrier=(17.6785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226574,'amu*angstrom^2'), symmetry=1, barrier=(17.6753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227189,'amu*angstrom^2'), symmetry=1, barrier=(17.6687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884634,0.0527144,4.89187e-06,-4.80413e-08,2.32242e-11,-310.197,27.4921], Tmin=(100,'K'), Tmax=(966.665,'K')), NASAPolynomial(coeffs=[17.3184,0.0227844,-7.74229e-06,1.415e-09,-1.03277e-13,-5266.19,-60.435], Tmin=(966.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.63004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]CCC=C=O(13151)',
    structure = SMILES('C[CH]CCC=C=O'),
    E0 = (45.2024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2120,512.5,787.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,304.146,304.151,1733.53],'cm^-1')),
        HinderedRotor(inertia=(0.00182263,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0978768,'amu*angstrom^2'), symmetry=1, barrier=(6.42415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0978613,'amu*angstrom^2'), symmetry=1, barrier=(6.42414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0978679,'amu*angstrom^2'), symmetry=1, barrier=(6.42411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55487,0.0583349,-4.32032e-05,2.04937e-08,-4.70718e-12,5520.78,26.8304], Tmin=(100,'K'), Tmax=(944.113,'K')), NASAPolynomial(coeffs=[5.08934,0.0433604,-1.94121e-05,3.69433e-09,-2.58796e-13,4853.38,9.98184], Tmin=(944.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.2024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(RCCJC)"""),
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
    label = 'CC[CH]CC=C[O](979)',
    structure = SMILES('CC[CH]CC=C[O]'),
    E0 = (63.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.790296,'amu*angstrom^2'), symmetry=1, barrier=(18.1705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00607604,'amu*angstrom^2'), symmetry=1, barrier=(3.63818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00279229,'amu*angstrom^2'), symmetry=1, barrier=(18.1862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0305484,'amu*angstrom^2'), symmetry=1, barrier=(18.1174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777747,0.0615351,-3.02277e-05,-9.00084e-10,3.8374e-12,7706.61,29.9689], Tmin=(100,'K'), Tmax=(1078.44,'K')), NASAPolynomial(coeffs=[12.9805,0.032601,-1.26922e-05,2.29791e-09,-1.58223e-13,4125.22,-34.2254], Tmin=(1078.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]CCCC=C[O](13152)',
    structure = SMILES('[CH2]CCCC=C[O]'),
    E0 = (73.8341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,430.845,430.908,430.922,430.949],'cm^-1')),
        HinderedRotor(inertia=(0.000907524,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113916,'amu*angstrom^2'), symmetry=1, barrier=(15.0121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113906,'amu*angstrom^2'), symmetry=1, barrier=(15.0122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113913,'amu*angstrom^2'), symmetry=1, barrier=(15.0122,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598614,0.0614135,-1.57489e-05,-2.54768e-08,1.46954e-11,9014.55,29.8313], Tmin=(100,'K'), Tmax=(986.983,'K')), NASAPolynomial(coeffs=[16.5713,0.0271154,-9.87854e-06,1.8016e-09,-1.28016e-13,4379.18,-54.5284], Tmin=(986.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.8341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CCC=[C]O(13153)',
    structure = SMILES('C[CH]CCC=[C]O'),
    E0 = (161.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3010,987.5,1337.5,450,1655,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519953,0.0700289,-5.62411e-05,2.45709e-08,-4.39539e-12,19532.4,32.7831], Tmin=(100,'K'), Tmax=(1325.75,'K')), NASAPolynomial(coeffs=[13.6999,0.0302628,-1.12482e-05,1.94568e-09,-1.28887e-13,16037.8,-34.5189], Tmin=(1325.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(RCCJC)"""),
)

species(
    label = 'CCC[CH]C=C[O](13154)',
    structure = SMILES('CCC[CH]C=C[O]'),
    E0 = (9.70036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767941,0.0553576,3.98887e-06,-4.66279e-08,2.22692e-11,1297.12,26.3304], Tmin=(100,'K'), Tmax=(976.653,'K')), NASAPolynomial(coeffs=[16.8126,0.0269376,-9.63845e-06,1.77116e-09,-1.27713e-13,-3615.47,-59.8022], Tmin=(976.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.70036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]CC[C]=CO(13155)',
    structure = SMILES('C[CH]CC[C]=CO'),
    E0 = (159.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3010,987.5,1337.5,450,1655,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331868,0.0722896,-5.24917e-05,1.35865e-08,1.14259e-12,19312.4,30.8], Tmin=(100,'K'), Tmax=(991.049,'K')), NASAPolynomial(coeffs=[15.5653,0.0273988,-9.66177e-06,1.66954e-09,-1.12943e-13,15478.2,-46.6665], Tmin=(991.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(RCCJC)"""),
)

species(
    label = 'CCCC[C]=C[O](13156)',
    structure = SMILES('CCCC[C]=C[O]'),
    E0 = (106.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,325.592,325.651,325.666,325.671],'cm^-1')),
        HinderedRotor(inertia=(0.00158936,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196028,'amu*angstrom^2'), symmetry=1, barrier=(14.7389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195976,'amu*angstrom^2'), symmetry=1, barrier=(14.7391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196009,'amu*angstrom^2'), symmetry=1, barrier=(14.7387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510082,0.0654886,-3.07766e-05,-7.6326e-09,7.84914e-12,12936,28.8366], Tmin=(100,'K'), Tmax=(1012.66,'K')), NASAPolynomial(coeffs=[15.7455,0.0286471,-1.07752e-05,1.95828e-09,-1.37077e-13,8653.69,-50.7658], Tmin=(1012.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C[CH]C[CH]C=CO(13157)',
    structure = SMILES('C[CH]C[CH]C=CO'),
    E0 = (62.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.590891,0.0621313,-1.75587e-05,-2.57444e-08,1.57645e-11,7673.5,28.2905], Tmin=(100,'K'), Tmax=(957.282,'K')), NASAPolynomial(coeffs=[16.6912,0.0255944,-8.47236e-06,1.47033e-09,-1.02599e-13,3182.6,-56.0373], Tmin=(957.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJC) + radical(Allyl_S)"""),
)

species(
    label = 'CCCC[CH][C]=O(4478)',
    structure = SMILES('CCCC[CH][C]=O'),
    E0 = (52.8261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89639,0.0657001,-4.75787e-05,1.88263e-08,-3.1477e-12,6467.12,27.8423], Tmin=(100,'K'), Tmax=(1367.17,'K')), NASAPolynomial(coeffs=[11.2897,0.0352921,-1.42165e-05,2.55813e-09,-1.7293e-13,3625.24,-25.5497], Tmin=(1367.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.8261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = 'C[CH][CH]CC=CO(13158)',
    structure = SMILES('C[CH][CH]CC=CO'),
    E0 = (116.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.524574,0.0691948,-5.48015e-05,2.37789e-08,-4.22753e-12,14086.3,32.2026], Tmin=(100,'K'), Tmax=(1335.01,'K')), NASAPolynomial(coeffs=[13.5685,0.0301127,-1.08898e-05,1.85083e-09,-1.21222e-13,10603.5,-34.4956], Tmin=(1335.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CCC=CO(13159)',
    structure = SMILES('[CH2][CH]CCC=CO'),
    E0 = (126.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424053,0.0681609,-3.72225e-05,-4.65763e-09,8.20171e-12,15390.8,31.7822], Tmin=(100,'K'), Tmax=(964.056,'K')), NASAPolynomial(coeffs=[16.4224,0.0258177,-8.73824e-06,1.50679e-09,-1.03395e-13,11189.2,-50.6081], Tmin=(964.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]C(44)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (279.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.00418548,'amu*angstrom^2'), symmetry=1, barrier=(6.91848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00418537,'amu*angstrom^2'), symmetry=1, barrier=(6.91838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25505,0.0137285,1.00536e-05,-1.43788e-08,4.3875e-12,33590.4,14.1736], Tmin=(100,'K'), Tmax=(1201.86,'K')), NASAPolynomial(coeffs=[3.74312,0.0203097,-8.40105e-06,1.5386e-09,-1.05137e-13,32880.4,9.26373], Tmin=(1201.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJC)"""),
)

species(
    label = 'C[CH][CH]CC=C[O](13160)',
    structure = SMILES('C[CH][CH]CC=C[O]'),
    E0 = (257.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,345.937,346.212,346.269,4000],'cm^-1')),
        HinderedRotor(inertia=(0.132132,'amu*angstrom^2'), symmetry=1, barrier=(11.2573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0514974,'amu*angstrom^2'), symmetry=1, barrier=(113.718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364209,'amu*angstrom^2'), symmetry=1, barrier=(31.0072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100228,'amu*angstrom^2'), symmetry=1, barrier=(113.799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24478,0.0586501,-3.90956e-05,1.41396e-08,-2.18297e-12,31069.6,31.0224], Tmin=(100,'K'), Tmax=(1447.21,'K')), NASAPolynomial(coeffs=[10.0951,0.0341883,-1.37415e-05,2.46002e-09,-1.65359e-13,28508,-14.9466], Tmin=(1447.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJC) + radical(C=COJ) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C[CH]C(38)',
    structure = SMILES('[CH2]C[CH]C'),
    E0 = (255.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1553.58],'cm^-1')),
        HinderedRotor(inertia=(0.00260972,'amu*angstrom^2'), symmetry=1, barrier=(4.49721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195012,'amu*angstrom^2'), symmetry=1, barrier=(4.4837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194147,'amu*angstrom^2'), symmetry=1, barrier=(4.46382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59222,0.0286125,-6.18646e-06,-3.09341e-09,1.21712e-12,30768.7,19.1945], Tmin=(100,'K'), Tmax=(1492.37,'K')), NASAPolynomial(coeffs=[6.22735,0.0257399,-1.0205e-05,1.78669e-09,-1.17174e-13,28918.6,-2.36167], Tmin=(1492.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC)"""),
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
    label = 'C[CH]C[CH]C=C[O](13161)',
    structure = SMILES('C[CH]C[CH]C=C[O]'),
    E0 = (204.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,679.847,679.898],'cm^-1')),
        HinderedRotor(inertia=(0.0883629,'amu*angstrom^2'), symmetry=1, barrier=(2.03164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00619457,'amu*angstrom^2'), symmetry=1, barrier=(2.03215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0393894,'amu*angstrom^2'), symmetry=1, barrier=(12.9192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118868,'amu*angstrom^2'), symmetry=1, barrier=(38.9886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986171,0.05533,-1.46224e-05,-1.91373e-08,1.09647e-11,24671.3,28.283], Tmin=(100,'K'), Tmax=(1006.31,'K')), NASAPolynomial(coeffs=[13.6213,0.0289778,-1.0924e-05,1.98509e-09,-1.38914e-13,20919.6,-38.7591], Tmin=(1006.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(C=COJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CCC=C[O](12802)',
    structure = SMILES('[CH2][CH]CCC=C[O]'),
    E0 = (268.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180,856.575,856.608,4000],'cm^-1')),
        HinderedRotor(inertia=(0.173442,'amu*angstrom^2'), symmetry=1, barrier=(3.98777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821541,'amu*angstrom^2'), symmetry=1, barrier=(18.8889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821569,'amu*angstrom^2'), symmetry=1, barrier=(18.8895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208278,'amu*angstrom^2'), symmetry=1, barrier=(108.485,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793707,0.0616451,-3.51915e-05,2.9619e-09,3.05205e-12,32389.7,31.8677], Tmin=(100,'K'), Tmax=(1050.09,'K')), NASAPolynomial(coeffs=[13.5649,0.0288505,-1.0992e-05,1.97555e-09,-1.3594e-13,28833.5,-34.5322], Tmin=(1050.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]CC[C]=C[O](13162)',
    structure = SMILES('C[CH]CC[C]=C[O]'),
    E0 = (300.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,180,393.645,393.676,393.676],'cm^-1')),
        HinderedRotor(inertia=(0.0199316,'amu*angstrom^2'), symmetry=1, barrier=(2.19208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13307,'amu*angstrom^2'), symmetry=1, barrier=(14.6386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13312,'amu*angstrom^2'), symmetry=1, barrier=(14.6386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00358897,'amu*angstrom^2'), symmetry=1, barrier=(40.7494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684811,0.0658934,-5.04683e-05,2.05536e-08,-3.42201e-12,36312.1,30.9499], Tmin=(100,'K'), Tmax=(1414.6,'K')), NASAPolynomial(coeffs=[13.8792,0.0285845,-1.09072e-05,1.90959e-09,-1.27094e-13,32579.1,-37.2818], Tmin=(1414.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(RCCJC) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]CC[CH][C]=O(13163)',
    structure = SMILES('C[CH]CC[CH][C]=O'),
    E0 = (247.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,245.042,884.453,1960.33],'cm^-1')),
        HinderedRotor(inertia=(0.100114,'amu*angstrom^2'), symmetry=1, barrier=(3.45534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100114,'amu*angstrom^2'), symmetry=1, barrier=(3.45534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100114,'amu*angstrom^2'), symmetry=1, barrier=(3.45534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100114,'amu*angstrom^2'), symmetry=1, barrier=(3.45534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100114,'amu*angstrom^2'), symmetry=1, barrier=(3.45534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25702,0.0638376,-5.92344e-05,3.69052e-08,-1.04096e-11,29835.6,29.2967], Tmin=(100,'K'), Tmax=(825.947,'K')), NASAPolynomial(coeffs=[5.99824,0.0408775,-1.75387e-05,3.25214e-09,-2.23901e-13,29052.4,7.32945], Tmin=(825.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(RCCJC) + radical(CCJCHO)"""),
)

species(
    label = 'C[CH]CCC1[CH]O1(9439)',
    structure = SMILES('C[CH]CCC1[CH]O1'),
    E0 = (203.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.197518,0.0711994,-5.79275e-05,2.70652e-08,-5.0608e-12,24610.5,30.0769], Tmin=(100,'K'), Tmax=(1445.24,'K')), NASAPolynomial(coeffs=[13.2236,0.0290607,-7.87508e-06,1.0629e-09,-5.88098e-14,21481,-35.3643], Tmin=(1445.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJC) + radical(CCsJO)"""),
)

species(
    label = 'CC1CC[CH]C1[O](13092)',
    structure = SMILES('CC1CC[CH]C1[O]'),
    E0 = (133.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84551,0.0233452,9.20989e-05,-1.35757e-07,5.35604e-11,16131.3,25.2752], Tmin=(100,'K'), Tmax=(951.684,'K')), NASAPolynomial(coeffs=[15.6441,0.0261098,-8.02772e-06,1.47073e-09,-1.11571e-13,10753.4,-55.0683], Tmin=(951.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC=CCC=CO(13164)',
    structure = SMILES('CC=CCC=CO'),
    E0 = (-156.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526984,0.0597225,-2.43736e-06,-4.76797e-08,2.50182e-11,-18630.5,26.9539], Tmin=(100,'K'), Tmax=(945.683,'K')), NASAPolynomial(coeffs=[19.5472,0.0206061,-5.95534e-06,1.0193e-09,-7.41938e-14,-24076.2,-73.5177], Tmin=(945.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'CCCCC=C=O(4486)',
    structure = SMILES('CCCCC=C=O'),
    E0 = (-149.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.977141,0.0627172,-4.00801e-05,1.29895e-08,-1.73961e-12,-17838.5,26.156], Tmin=(100,'K'), Tmax=(1671.91,'K')), NASAPolynomial(coeffs=[13.2885,0.0332628,-1.36546e-05,2.45255e-09,-1.64045e-13,-21955.2,-39.5671], Tmin=(1671.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'C=CCCC=CO(13165)',
    structure = SMILES('C=CCCC=CO'),
    E0 = (-145.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.474325,0.0596763,1.48864e-06,-5.42134e-08,2.79288e-11,-17307.3,27.5544], Tmin=(100,'K'), Tmax=(941.884,'K')), NASAPolynomial(coeffs=[20.5442,0.0191422,-5.1439e-06,8.66025e-10,-6.4151e-14,-23070.7,-78.5947], Tmin=(941.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[CH]C[CH]C[CH][O](13166)',
    structure = SMILES('C[CH]C[CH]C[CH][O]'),
    E0 = (443.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,314.019,792.813,878.806,2884.31,3332.8],'cm^-1')),
        HinderedRotor(inertia=(0.0581298,'amu*angstrom^2'), symmetry=1, barrier=(1.33652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0581298,'amu*angstrom^2'), symmetry=1, barrier=(1.33652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0581298,'amu*angstrom^2'), symmetry=1, barrier=(1.33652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0581298,'amu*angstrom^2'), symmetry=1, barrier=(1.33652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0581298,'amu*angstrom^2'), symmetry=1, barrier=(1.33652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.91445,0.0836024,-0.000137704,1.4232e-07,-5.53372e-11,53402.5,33.5784], Tmin=(100,'K'), Tmax=(859.786,'K')), NASAPolynomial(coeffs=[-2.11768,0.0592392,-2.80848e-05,5.28257e-09,-3.59111e-13,55345.8,56.0175], Tmin=(859.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(RCCJCC) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH][CH]CC[CH][O](992)',
    structure = SMILES('C[CH][CH]CC[CH][O]'),
    E0 = (443.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,180,849.253,1000.5,2720.92,3515.24],'cm^-1')),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.914328,0.083604,-0.00013771,1.42329e-07,-5.53415e-11,53402.5,33.5788], Tmin=(100,'K'), Tmax=(859.77,'K')), NASAPolynomial(coeffs=[-2.11741,0.0592388,-2.80846e-05,5.2825e-09,-3.59106e-13,55345.7,56.016], Tmin=(859.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH][CH]C[CH]C[O](13167)',
    structure = SMILES('C[CH][CH]C[CH]C[O]'),
    E0 = (462.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,512.234,695.237,2737.89,3355.56,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0625048,'amu*angstrom^2'), symmetry=1, barrier=(2.49573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0625048,'amu*angstrom^2'), symmetry=1, barrier=(2.49573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0625048,'amu*angstrom^2'), symmetry=1, barrier=(2.49573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0625048,'amu*angstrom^2'), symmetry=1, barrier=(2.49573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0625048,'amu*angstrom^2'), symmetry=1, barrier=(2.49573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35062,0.0699868,-9.91344e-05,1.02359e-07,-4.11967e-11,55748.4,34.9598], Tmin=(100,'K'), Tmax=(840.954,'K')), NASAPolynomial(coeffs=[-2.10693,0.0580293,-2.71434e-05,5.12514e-09,-3.51255e-13,57334.3,57.013], Tmin=(840.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCOJ) + radical(CCJCO) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C(C)CC=C[O](12592)',
    structure = SMILES('[CH2]C(C)CC=C[O]'),
    E0 = (64.6285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,830.114,831.194],'cm^-1')),
        HinderedRotor(inertia=(0.146011,'amu*angstrom^2'), symmetry=1, barrier=(15.8569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.68954,'amu*angstrom^2'), symmetry=1, barrier=(15.8539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69015,'amu*angstrom^2'), symmetry=1, barrier=(15.8679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0324472,'amu*angstrom^2'), symmetry=1, barrier=(15.8655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3875.05,'J/mol'), sigma=(6.66255,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.27 K, Pc=29.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649327,0.0591604,-5.27715e-06,-4.08907e-08,2.18549e-11,7906.92,29.3645], Tmin=(100,'K'), Tmax=(942.813,'K')), NASAPolynomial(coeffs=[17.3392,0.0244432,-7.46431e-06,1.25881e-09,-8.80823e-14,3155.72,-58.6788], Tmin=(942.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.6285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = 'CC1CCC=CO1(13034)',
    structure = SMILES('CC1CCC=CO1'),
    E0 = (-190.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56378,0.0268899,9.54199e-05,-1.50664e-07,6.24109e-11,-22836.7,18.1438], Tmin=(100,'K'), Tmax=(925.397,'K')), NASAPolynomial(coeffs=[19.2107,0.0191688,-3.19117e-06,4.33233e-10,-3.6262e-14,-29038.3,-81.4849], Tmin=(925.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran)"""),
)

species(
    label = '[CH2]CC=C[O](743)',
    structure = SMILES('[CH2]CC=C[O]'),
    E0 = (121.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,439.352,440.116],'cm^-1')),
        HinderedRotor(inertia=(0.132298,'amu*angstrom^2'), symmetry=1, barrier=(18.3159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132458,'amu*angstrom^2'), symmetry=1, barrier=(18.2855,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06392,0.0297838,1.93326e-05,-5.14894e-08,2.33454e-11,14681.8,20.083], Tmin=(100,'K'), Tmax=(942.939,'K')), NASAPolynomial(coeffs=[13.9169,0.0116798,-3.05432e-06,5.27558e-10,-4.05895e-14,11016,-43.9895], Tmin=(942.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]C(32)',
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
    label = '[CH3](11)',
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
    label = '[CH]CCC=C[O](12968)',
    structure = SMILES('[CH]CCC=C[O]'),
    E0 = (340.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,497.853,497.872,497.875,497.9,497.91,497.929],'cm^-1')),
        HinderedRotor(inertia=(0.000680023,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086622,'amu*angstrom^2'), symmetry=1, barrier=(15.238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866422,'amu*angstrom^2'), symmetry=1, barrier=(15.2374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25001,0.0478266,-6.97292e-06,-3.13923e-08,1.71937e-11,41073.2,24.4661], Tmin=(100,'K'), Tmax=(955.545,'K')), NASAPolynomial(coeffs=[16.2596,0.01576,-4.92966e-06,8.76369e-10,-6.42451e-14,36800.2,-54.6134], Tmin=(955.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCJ2_triplet) + radical(C=COJ)"""),
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
    label = '[CH]=CCC[CH]C(5241)',
    structure = SMILES('[CH]=CCC[CH]C'),
    E0 = (377.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,346.469,1300.99],'cm^-1')),
        HinderedRotor(inertia=(0.0702507,'amu*angstrom^2'), symmetry=1, barrier=(5.98546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702519,'amu*angstrom^2'), symmetry=1, barrier=(5.98548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702466,'amu*angstrom^2'), symmetry=1, barrier=(5.98549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702462,'amu*angstrom^2'), symmetry=1, barrier=(5.98552,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58457,0.0524453,-2.90832e-05,8.01022e-09,-9.1372e-13,45484.6,26.3683], Tmin=(100,'K'), Tmax=(1875.64,'K')), NASAPolynomial(coeffs=[11.0482,0.0322633,-1.29432e-05,2.27355e-09,-1.49096e-13,41934.5,-25.2403], Tmin=(1875.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_P)"""),
)

species(
    label = 'C[C]CCC=C[O](13168)',
    structure = SMILES('C[C]CCC=C[O]'),
    E0 = (316.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.514118,0.0636794,-2.46438e-05,-1.82243e-08,1.28044e-11,38239.7,28.2532], Tmin=(100,'K'), Tmax=(979.267,'K')), NASAPolynomial(coeffs=[17.6047,0.0234474,-8.32434e-06,1.50928e-09,-1.07618e-13,33474.2,-61.0821], Tmin=(979.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH]CC=CC=O(13169)',
    structure = SMILES('C[CH]CC=CC=O'),
    E0 = (30.6562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,246.588,246.589],'cm^-1')),
        HinderedRotor(inertia=(0.00277238,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277223,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160495,'amu*angstrom^2'), symmetry=1, barrier=(6.92555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160502,'amu*angstrom^2'), symmetry=1, barrier=(6.92556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09505,0.0523779,-2.73552e-05,6.11035e-09,-5.16978e-13,3744.25,24.9236], Tmin=(100,'K'), Tmax=(2729.76,'K')), NASAPolynomial(coeffs=[25.6559,0.0178534,-8.38406e-06,1.4772e-09,-9.26596e-14,-9118.88,-112.404], Tmin=(2729.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.6562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C[CH]CC=O(13120)',
    structure = SMILES('C[CH]C[CH]CC=O'),
    E0 = (119.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,214.751,398.094,3427.35],'cm^-1')),
        HinderedRotor(inertia=(0.0396268,'amu*angstrom^2'), symmetry=1, barrier=(1.22942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0396268,'amu*angstrom^2'), symmetry=1, barrier=(1.22942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0396268,'amu*angstrom^2'), symmetry=1, barrier=(1.22942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0396268,'amu*angstrom^2'), symmetry=1, barrier=(1.22942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0396268,'amu*angstrom^2'), symmetry=1, barrier=(1.22942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22834,0.047793,2.49037e-05,-1.37136e-07,1.22105e-10,14449.2,27.9365], Tmin=(100,'K'), Tmax=(427.851,'K')), NASAPolynomial(coeffs=[3.46105,0.0487361,-2.21137e-05,4.2348e-09,-2.97787e-13,14229.6,21.7023], Tmin=(427.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]CCC[C]=O(4477)',
    structure = SMILES('C[CH]CCC[C]=O'),
    E0 = (79.743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01064,0.0713114,-8.13435e-05,6.67281e-08,-2.38944e-11,9693.16,30.411], Tmin=(100,'K'), Tmax=(790.948,'K')), NASAPolynomial(coeffs=[4.33771,0.0470036,-2.10552e-05,3.95303e-09,-2.72468e-13,9400.9,16.6195], Tmin=(790.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(RCCJC)"""),
)

species(
    label = 'CC=CCC[CH][O](978)',
    structure = SMILES('CC=CCC[CH][O]'),
    E0 = (169.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,259.53,259.534,259.563,2103.12],'cm^-1')),
        HinderedRotor(inertia=(0.00249932,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179356,'amu*angstrom^2'), symmetry=1, barrier=(8.57114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179246,'amu*angstrom^2'), symmetry=1, barrier=(8.56682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17928,'amu*angstrom^2'), symmetry=1, barrier=(8.5686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.88854,0.0744443,-8.59456e-05,7.09703e-08,-2.58986e-11,20534.3,28.4129], Tmin=(100,'K'), Tmax=(757.268,'K')), NASAPolynomial(coeffs=[4.57547,0.0487226,-2.26224e-05,4.3301e-09,-3.02226e-13,20155,12.8333], Tmin=(757.268,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]CCCC=O(13170)',
    structure = SMILES('[CH2][CH]CCCC=O'),
    E0 = (125.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,190.234,322.697,3062.46],'cm^-1')),
        HinderedRotor(inertia=(0.064475,'amu*angstrom^2'), symmetry=1, barrier=(1.62689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.064475,'amu*angstrom^2'), symmetry=1, barrier=(1.62689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.064475,'amu*angstrom^2'), symmetry=1, barrier=(1.62689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.064475,'amu*angstrom^2'), symmetry=1, barrier=(1.62689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.064475,'amu*angstrom^2'), symmetry=1, barrier=(1.62689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13378,0.0689471,-7.64494e-05,6.47162e-08,-2.43681e-11,15135,30.936], Tmin=(100,'K'), Tmax=(775.668,'K')), NASAPolynomial(coeffs=[3.17539,0.0498192,-2.28295e-05,4.33809e-09,-3.01329e-13,15077,23.2726], Tmin=(775.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'CC1CC[CH][CH]O1(13041)',
    structure = SMILES('CC1CC[CH][CH]O1'),
    E0 = (99.7724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4621,0.035409,6.55354e-05,-1.17375e-07,5.10493e-11,12110.4,21.7666], Tmin=(100,'K'), Tmax=(901.3,'K')), NASAPolynomial(coeffs=[16.3257,0.0233308,-4.04583e-06,4.27766e-10,-2.65803e-14,7242.29,-60.5394], Tmin=(901.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.7724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxane) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = 'CC=CCCC=O(4487)',
    structure = SMILES('CC=CCCC=O'),
    E0 = (-159.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4883,0.0578452,-3.2577e-05,8.7911e-09,-9.69802e-13,-19051.8,25.1553], Tmin=(100,'K'), Tmax=(1926.54,'K')), NASAPolynomial(coeffs=[12.3753,0.0352411,-1.49777e-05,2.70104e-09,-1.79526e-13,-23246.7,-34.5074], Tmin=(1926.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCCC=CC=O(13171)',
    structure = SMILES('CCCC=CC=O'),
    E0 = (-163.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04072,0.0612498,-3.58396e-05,9.56138e-09,-9.9796e-13,-19589.9,26.0521], Tmin=(100,'K'), Tmax=(2186.04,'K')), NASAPolynomial(coeffs=[20.9614,0.0247985,-1.08275e-05,1.93344e-09,-1.25603e-13,-28299.3,-85.6335], Tmin=(2186.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=CCCCC=O(13172)',
    structure = SMILES('C=CCCCC=O'),
    E0 = (-146.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21803,0.060042,-3.55964e-05,1.03372e-08,-1.23118e-12,-17565.6,26.5806], Tmin=(100,'K'), Tmax=(1834.9,'K')), NASAPolynomial(coeffs=[13.2573,0.0337967,-1.41412e-05,2.54192e-09,-1.69092e-13,-21983.8,-38.8097], Tmin=(1834.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(C=O)C[CH]C(12590)',
    structure = SMILES('[CH2]C(C=O)C[CH]C'),
    E0 = (124.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,180,1619.74],'cm^-1')),
        HinderedRotor(inertia=(0.177767,'amu*angstrom^2'), symmetry=1, barrier=(4.08721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17787,'amu*angstrom^2'), symmetry=1, barrier=(4.08958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177882,'amu*angstrom^2'), symmetry=1, barrier=(4.08986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177849,'amu*angstrom^2'), symmetry=1, barrier=(4.0891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1779,'amu*angstrom^2'), symmetry=1, barrier=(4.09026,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80362,0.07634,-9.20439e-05,7.62615e-08,-2.70757e-11,15087.1,30.5484], Tmin=(100,'K'), Tmax=(792.621,'K')), NASAPolynomial(coeffs=[5.05372,0.0469085,-2.12382e-05,4.00038e-09,-2.75984e-13,14664.1,12.6137], Tmin=(792.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CC1CCC1C=O(12596)',
    structure = SMILES('CC1CCC1C=O'),
    E0 = (-131.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60327,0.0358881,4.53963e-05,-7.82496e-08,3.05707e-11,-15689.2,23.1854], Tmin=(100,'K'), Tmax=(994.35,'K')), NASAPolynomial(coeffs=[12.4534,0.0333815,-1.28835e-05,2.43351e-09,-1.76086e-13,-19880.8,-39.3256], Tmin=(994.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH]CC[CH]C(131)',
    structure = SMILES('[CH]CC[CH]C'),
    E0 = (474.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,206.934,850.427,1133.9,1446.07,1744.2],'cm^-1')),
        HinderedRotor(inertia=(0.131703,'amu*angstrom^2'), symmetry=1, barrier=(3.2876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131703,'amu*angstrom^2'), symmetry=1, barrier=(3.2876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131703,'amu*angstrom^2'), symmetry=1, barrier=(3.2876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131703,'amu*angstrom^2'), symmetry=1, barrier=(3.2876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08309,0.0434474,-2.26958e-05,5.77488e-09,-6.05056e-13,57145.3,22.4528], Tmin=(100,'K'), Tmax=(1992.05,'K')), NASAPolynomial(coeffs=[9.58095,0.0283916,-1.13587e-05,1.9807e-09,-1.28884e-13,54158.2,-18.8872], Tmin=(1992.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCJ2_triplet)"""),
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
    E0 = (63.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (194.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (202.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (210.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (259.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (142.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (222.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (221.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (324.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (190.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (309.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (150.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (176.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (89.9162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (210.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (207.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (369.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (469.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (478.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (421.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (480.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (513.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (459.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (276.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (133.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (97.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (71.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (88.0072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (465.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (506.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (487.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (224.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (70.5652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (507.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (510.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (784.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (528.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (242.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (202.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (244.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (239.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (294.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (184.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (99.7724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (162.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (126.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (102.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (284.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (71.3183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (542.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['C=CC(42)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CC1CCC1[CH][O](13149)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(966131,'s^-1'), n=1.86605, Ea=(131.371,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 129.8 to 131.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'CC=CCC=C[O](13150)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C=CCCC=C[O](12800)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C[CH]CCC=C=O(13151)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=CC(42)', '[CH2]C=C[O](5266)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(780,'cm^3/(mol*s)','*|/',2), n=2.53, Ea=(46.024,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 4 used for Cds-HH_Cds-Cs\H3/H;CsJ-CdHH
Exact match found for rate rule [Cds-HH_Cds-Cs\H3/H;CsJ-CdHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC[CH]CC=C[O](979)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]CCCC=C[O](13152)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH]CCC=[C]O(13153)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CCC[CH]C=C[O](13154)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 219 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]CC[C]=CO(13155)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCCC[C]=C[O](13156)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['C[CH]C[CH]C=CO(13157)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 289 used for R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CCCC[CH][C]=O(4478)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSR;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_SSSD;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH][CH]CC=CO(13158)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]CCC=CO(13159)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(60975.7,'s^-1'), n=1.58648, Ea=(80.9836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_2H;XH_out] for rate rule [R7HJ_1;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C(44)', '[CH2]C=C[O](5266)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', 'C[CH][CH]CC=C[O](13160)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C[CH]C(38)', '[CH]=C[O](602)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.00218e+08,'m^3/(mol*s)'), n=-0.446058, Ea=(0.74957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H2/Cs] + [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/Cs]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', 'C[CH]C[CH]C=C[O](13161)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.66284e+07,'m^3/(mol*s)'), n=0.108445, Ea=(5.74998,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2][CH]CCC=C[O](12802)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', 'C[CH]CC[C]=C[O](13162)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', 'C[CH]CC[CH][C]=O(13163)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['C[CH]CCC1[CH]O1(9439)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CC1CC[CH]C1[O](13092)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(70.2628,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_CsCs_HH_D;doublebond_intra_pri;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 65.5 to 70.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CC=CCC=CO(13164)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CCCCC=C=O(4486)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['C=CCCC=CO(13165)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[CH]C[CH]C[CH][O](13166)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C[CH][CH]CC[CH][O](992)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH][CH]C[CH]C[O](13167)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C)CC=C[O](12592)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CC1CCC=CO1(13034)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_1H;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]CC=C[O](743)', '[CH]C(32)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH3](11)', '[CH]CCC=C[O](12968)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(T)(63)', '[CH]=CCC[CH]C(5241)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(8)', 'C[C]CCC=C[O](13168)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(8)', 'C[CH]CC=CC=O(13169)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.76955,'m^3/(mol*s)'), n=1.94497, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-OneDeH;HJ] for rate rule [Cds-CsH_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]C(44)', 'C=CC=O(5269)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.0102751,'m^3/(mol*s)'), n=2.40501, Ea=(4.48561,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDeH;CJ] for rate rule [Cds-HH_Cds-COH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C[CH]C[CH]CC=O(13120)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C[CH]CCC[C]=O(4477)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CC=CCC[CH][O](978)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][CH]CCCC=O(13170)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(61977.9,'s^-1'), n=1.86063, Ea=(59.8007,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CC1CC[CH][CH]O1(13041)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(36.7384,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCs] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 34.0 to 36.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction45',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CC=CCCC=O(4487)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.94203e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CCCC=CC=O(13171)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['C=CCCCC=O(13172)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C[CH]CCC=C[O](12588)'],
    products = ['CC1CCC1C=O(12596)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]CC[CH]C(131)', '[CH]=O(373)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '3322',
    isomers = [
        'C[CH]CCC=C[O](12588)',
    ],
    reactants = [
        ('C=CC(42)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3322',
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

