species(
    label = '[CH2][CH]OO[CH]C=C(12655)',
    structure = SMILES('[CH2][CH]OO[CH]C=C'),
    E0 = (404.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,350,500,795,815,3010,987.5,1337.5,450,1655,282.25,282.279],'cm^-1')),
        HinderedRotor(inertia=(0.535678,'amu*angstrom^2'), symmetry=1, barrier=(30.2929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0165703,'amu*angstrom^2'), symmetry=1, barrier=(30.2948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0267644,'amu*angstrom^2'), symmetry=1, barrier=(48.9303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.535942,'amu*angstrom^2'), symmetry=1, barrier=(30.2937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0882041,'amu*angstrom^2'), symmetry=1, barrier=(30.2937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702125,0.0652698,-5.16693e-05,2.0244e-08,-3.1752e-12,48782.8,30.5043], Tmin=(100,'K'), Tmax=(1506.08,'K')), NASAPolynomial(coeffs=[16.3679,0.0236638,-1.02321e-05,1.90216e-09,-1.30623e-13,44064,-51.4895], Tmin=(1506.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(CJCOOH) + radical(CCsJOOC)"""),
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
    label = '[CH2]C1[CH]OOC1[CH2](14687)',
    structure = SMILES('[CH2]C1[CH]OOC1[CH2]'),
    E0 = (411.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3572,0.0444299,1.8094e-05,-6.87111e-08,3.55215e-11,49638.9,24.9957], Tmin=(100,'K'), Tmax=(862.85,'K')), NASAPolynomial(coeffs=[16.4519,0.0147456,-3.46148e-07,-3.45243e-10,3.32496e-14,45534.1,-54.2924], Tmin=(862.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(Isobutyl) + radical(CCsJOOC) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C1[CH]OO[CH]C1(14688)',
    structure = SMILES('[CH2]C1[CH]OO[CH]C1'),
    E0 = (388.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59366,0.0389162,3.08012e-05,-7.76804e-08,3.74116e-11,46848.2,21.5888], Tmin=(100,'K'), Tmax=(870.263,'K')), NASAPolynomial(coeffs=[14.7654,0.0181619,-2.0039e-06,-1.59928e-11,9.57672e-15,43049,-48.7831], Tmin=(870.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(CCsJOOC) + radical(Isobutyl) + radical(CCsJOOC)"""),
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
    label = 'C=[C][CH]OOC=C(14336)',
    structure = SMILES('C=[C][CH]OOC=C'),
    E0 = (384.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,500,795,815,3010,987.5,1337.5,450,1655,405.171,405.172],'cm^-1')),
        HinderedRotor(inertia=(0.234103,'amu*angstrom^2'), symmetry=1, barrier=(27.2716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234103,'amu*angstrom^2'), symmetry=1, barrier=(27.2716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234103,'amu*angstrom^2'), symmetry=1, barrier=(27.2716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234103,'amu*angstrom^2'), symmetry=1, barrier=(27.2716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.962042,0.0559893,-2.92671e-05,-5.08293e-09,6.38531e-12,46414.4,28.4269], Tmin=(100,'K'), Tmax=(1035.58,'K')), NASAPolynomial(coeffs=[16.013,0.0196042,-8.06878e-06,1.55154e-09,-1.12394e-13,42130.8,-50.3425], Tmin=(1035.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][CH]OOC[C]=C(14689)',
    structure = SMILES('[CH2][CH]OOC[C]=C'),
    E0 = (525.114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,184.985,1681.58],'cm^-1')),
        HinderedRotor(inertia=(0.402819,'amu*angstrom^2'), symmetry=1, barrier=(9.78119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00487453,'amu*angstrom^2'), symmetry=1, barrier=(9.78118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.92635,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58635,'amu*angstrom^2'), symmetry=1, barrier=(38.5212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58635,'amu*angstrom^2'), symmetry=1, barrier=(38.5212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45698,0.0664728,-4.74143e-05,-3.43304e-08,5.85463e-11,63237.5,27.8169], Tmin=(100,'K'), Tmax=(488.63,'K')), NASAPolynomial(coeffs=[7.08113,0.0388029,-1.88662e-05,3.66042e-09,-2.56671e-13,62468.6,2.4673], Tmin=(488.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOOC) + radical(Cds_S) + radical(CJCOOH)"""),
)

species(
    label = '[CH]=CCOO[CH][CH2](14690)',
    structure = SMILES('[CH]=CCOO[CH][CH2]'),
    E0 = (534.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,269.466],'cm^-1')),
        HinderedRotor(inertia=(0.151586,'amu*angstrom^2'), symmetry=1, barrier=(7.81069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00345232,'amu*angstrom^2'), symmetry=1, barrier=(7.81069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754596,'amu*angstrom^2'), symmetry=1, barrier=(38.8806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754566,'amu*angstrom^2'), symmetry=1, barrier=(38.8806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.32164,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.881783,0.0743747,-8.92256e-05,6.40007e-08,-1.94863e-11,64376.8,29.6707], Tmin=(100,'K'), Tmax=(786.377,'K')), NASAPolynomial(coeffs=[8.23929,0.0369499,-1.78385e-05,3.48095e-09,-2.46214e-13,63219.6,-4.05688], Tmin=(786.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOOC) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]COO[CH][C]=C(14691)',
    structure = SMILES('[CH2]COO[CH][C]=C'),
    E0 = (455.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,250.255,667.786],'cm^-1')),
        HinderedRotor(inertia=(0.641849,'amu*angstrom^2'), symmetry=1, barrier=(28.5188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.640288,'amu*angstrom^2'), symmetry=1, barrier=(28.5175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(7.76333e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.640565,'amu*angstrom^2'), symmetry=1, barrier=(28.519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961981,'amu*angstrom^2'), symmetry=1, barrier=(42.7353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888183,0.0655049,-5.42026e-05,2.2689e-08,-3.85953e-12,54954.6,30.0049], Tmin=(100,'K'), Tmax=(1378.88,'K')), NASAPolynomial(coeffs=[14.2656,0.0266982,-1.19873e-05,2.27863e-09,-1.59008e-13,51265.4,-38.8314], Tmin=(1378.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCOOH) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]OO[CH]C(13626)',
    structure = SMILES('C=[C][CH]OO[CH]C'),
    E0 = (428.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,500,795,815,3000,3050,390,425,1340,1360,335,370,314.987,981.484],'cm^-1')),
        HinderedRotor(inertia=(0.473409,'amu*angstrom^2'), symmetry=1, barrier=(33.2171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471763,'amu*angstrom^2'), symmetry=1, barrier=(33.2175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0485883,'amu*angstrom^2'), symmetry=1, barrier=(33.2171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.719488,'amu*angstrom^2'), symmetry=1, barrier=(50.5565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0757392,'amu*angstrom^2'), symmetry=1, barrier=(5.33045,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11956,0.0628063,-4.89512e-05,1.93901e-08,-3.174e-12,51634,27.471], Tmin=(100,'K'), Tmax=(1402.3,'K')), NASAPolynomial(coeffs=[12.6888,0.0298054,-1.36511e-05,2.60805e-09,-1.82123e-13,48389.3,-32.2558], Tmin=(1402.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOOC) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]OOC[CH2](14692)',
    structure = SMILES('[CH]=C[CH]OOC[CH2]'),
    E0 = (465.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.37945,'amu*angstrom^2'), symmetry=1, barrier=(31.7162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0535238,'amu*angstrom^2'), symmetry=1, barrier=(31.7154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37941,'amu*angstrom^2'), symmetry=1, barrier=(31.7153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175911,'amu*angstrom^2'), symmetry=1, barrier=(4.04453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.053585,'amu*angstrom^2'), symmetry=1, barrier=(31.7204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576778,0.0674854,-5.61855e-05,2.31778e-08,-3.80642e-12,56083,31.0611], Tmin=(100,'K'), Tmax=(1451.98,'K')), NASAPolynomial(coeffs=[16.9481,0.0223846,-9.59297e-06,1.78515e-09,-1.23047e-13,51328.9,-54.0265], Tmin=(1451.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCOOH) + radical(Cds_P) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=C[CH]OO[CH]C(13627)',
    structure = SMILES('[CH]=C[CH]OO[CH]C'),
    E0 = (437.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,500,795,815,3000,3050,390,425,1340,1360,335,370,223.595],'cm^-1')),
        HinderedRotor(inertia=(0.0577519,'amu*angstrom^2'), symmetry=1, barrier=(29.8722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.832201,'amu*angstrom^2'), symmetry=1, barrier=(29.8742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104992,'amu*angstrom^2'), symmetry=1, barrier=(29.8758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842213,'amu*angstrom^2'), symmetry=1, barrier=(29.8766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50649,'amu*angstrom^2'), symmetry=1, barrier=(53.3928,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.806128,0.0648102,-5.10091e-05,1.99617e-08,-3.14928e-12,52762.5,28.5344], Tmin=(100,'K'), Tmax=(1486.57,'K')), NASAPolynomial(coeffs=[15.4795,0.0253277,-1.11697e-05,2.09527e-09,-1.44637e-13,48400,-48.0734], Tmin=(1486.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(Cds_P) + radical(CCsJOOC)"""),
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
    label = '[CH2][CH]OO[CH][C]=C(14693)',
    structure = SMILES('[CH2][CH]OO[CH][C]=C'),
    E0 = (642.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,350,500,795,815,3000,3100,440,815,1455,1000,300.234,1627.99],'cm^-1')),
        HinderedRotor(inertia=(0.0139285,'amu*angstrom^2'), symmetry=1, barrier=(26.2047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148336,'amu*angstrom^2'), symmetry=1, barrier=(52.7776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805962,'amu*angstrom^2'), symmetry=1, barrier=(52.7796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.498083,'amu*angstrom^2'), symmetry=1, barrier=(25.668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386054,'amu*angstrom^2'), symmetry=1, barrier=(25.6602,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11954,0.0647068,-5.91536e-05,2.82136e-08,-5.55471e-12,77366.5,29.4715], Tmin=(100,'K'), Tmax=(1192.5,'K')), NASAPolynomial(coeffs=[11.9803,0.0282767,-1.33297e-05,2.59576e-09,-1.84103e-13,74776.2,-24.8373], Tmin=(1192.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(642.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJO) + radical(CJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = '[CH]=C[CH]OO[CH][CH2](14694)',
    structure = SMILES('[CH]=C[CH]OO[CH][CH2]'),
    E0 = (651.664,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3000,3050,390,425,1340,1360,335,370,350,500,795,815,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,604.208],'cm^-1')),
        HinderedRotor(inertia=(0.0196586,'amu*angstrom^2'), symmetry=1, barrier=(28.4021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025375,'amu*angstrom^2'), symmetry=1, barrier=(36.5782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.72197,'amu*angstrom^2'), symmetry=1, barrier=(85.5753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23419,'amu*angstrom^2'), symmetry=1, barrier=(28.3764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23423,'amu*angstrom^2'), symmetry=1, barrier=(28.3774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.880591,0.0659709,-5.91784e-05,2.67761e-08,-4.88398e-12,78491.4,30.2572], Tmin=(100,'K'), Tmax=(1302.05,'K')), NASAPolynomial(coeffs=[14.4909,0.0241593,-1.10105e-05,2.11366e-09,-1.48699e-13,74947.2,-38.9968], Tmin=(1302.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(CJCOOH) + radical(Cds_P) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2][CH]OOC1[CH]C1(12044)',
    structure = SMILES('[CH2][CH]OOC1[CH]C1'),
    E0 = (504.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,350,500,795,815,3025,407.5,1350,352.5,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790149,0.0634706,-4.95753e-05,1.92705e-08,-3.00463e-12,60774.8,30.2465], Tmin=(100,'K'), Tmax=(1512.67,'K')), NASAPolynomial(coeffs=[15.8553,0.0236337,-1.00723e-05,1.86077e-09,-1.27324e-13,56217.1,-48.6692], Tmin=(1512.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CCJCOOH) + radical(CJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]C1C[CH][CH]OO1(14695)',
    structure = SMILES('[CH2]C1C[CH][CH]OO1'),
    E0 = (407.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43325,0.0414827,2.52934e-05,-7.20813e-08,3.47193e-11,49170.3,24.3192], Tmin=(100,'K'), Tmax=(890.889,'K')), NASAPolynomial(coeffs=[15.9837,0.0170388,-2.39064e-06,1.49515e-10,-5.87413e-15,44955.3,-53.3034], Tmin=(890.889,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(CJCOOH) + radical(CCJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = '[CH]1[CH]OO[CH]CC1(14696)',
    structure = SMILES('[CH]1[CH]OO[CH]CC1'),
    E0 = (406.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47415,0.0104101,0.000168143,-2.59048e-07,1.10222e-10,49083.1,26.1301], Tmin=(100,'K'), Tmax=(907.535,'K')), NASAPolynomial(coeffs=[33.9829,-0.0148098,1.46877e-05,-2.97402e-09,1.9268e-13,38320.5,-154.339], Tmin=(907.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(CCsJOOC) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]COOC=C=C(14697)',
    structure = SMILES('[CH2]COOC=C=C'),
    E0 = (274.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749666,0.0700401,-6.63569e-05,3.26908e-08,-6.51845e-12,33165.3,27.8706], Tmin=(100,'K'), Tmax=(1198.07,'K')), NASAPolynomial(coeffs=[13.8638,0.0262557,-1.1538e-05,2.18666e-09,-1.53141e-13,30023,-37.7673], Tmin=(1198.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C1OOC1C=C(12657)',
    structure = SMILES('[CH2]C1OOC1C=C'),
    E0 = (235.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3285,0.0414507,2.97277e-05,-8.10966e-08,3.89489e-11,28380.4,25.0467], Tmin=(100,'K'), Tmax=(895.02,'K')), NASAPolynomial(coeffs=[18.2883,0.0125156,-3.15714e-07,-2.19211e-10,1.77516e-14,23467.5,-65.379], Tmin=(895.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CJCOOH)"""),
)

species(
    label = 'C=C[CH]OOC=C(12635)',
    structure = SMILES('C=C[CH]OOC=C'),
    E0 = (147.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680,382.825,382.975],'cm^-1')),
        HinderedRotor(inertia=(0.268761,'amu*angstrom^2'), symmetry=1, barrier=(28.0167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269293,'amu*angstrom^2'), symmetry=1, barrier=(28.0253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268149,'amu*angstrom^2'), symmetry=1, barrier=(28.0214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268568,'amu*angstrom^2'), symmetry=1, barrier=(28.0146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01094,0.0512204,-3.98472e-06,-3.48326e-08,1.74765e-11,17810.4,27.7766], Tmin=(100,'K'), Tmax=(996.028,'K')), NASAPolynomial(coeffs=[17.0046,0.0204305,-7.97523e-06,1.54532e-09,-1.14572e-13,12965.6,-57.6465], Tmin=(996.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO)"""),
)

species(
    label = '[CH][CH2](721)',
    structure = SMILES('[CH][CH2]'),
    E0 = (556.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1101.59,1101.66],'cm^-1')),
        HinderedRotor(inertia=(0.00420677,'amu*angstrom^2'), symmetry=1, barrier=(3.62356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77493,-0.000462567,3.18167e-05,-4.30783e-08,1.77606e-11,66973.8,8.79001], Tmin=(100,'K'), Tmax=(870.354,'K')), NASAPolynomial(coeffs=[6.06996,0.00332438,5.85464e-07,-2.32999e-10,1.82455e-14,66031.4,-5.08252], Tmin=(870.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCJ)"""),
)

species(
    label = 'C=C[CH]O[O](6572)',
    structure = SMILES('C=C[CH]O[O]'),
    E0 = (193.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.0263592,'amu*angstrom^2'), symmetry=1, barrier=(16.4069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264447,'amu*angstrom^2'), symmetry=1, barrier=(37.5506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44754,0.0263949,2.69871e-06,-2.29929e-08,1.07611e-11,23370.7,18.6111], Tmin=(100,'K'), Tmax=(985.371,'K')), NASAPolynomial(coeffs=[10.196,0.0131636,-4.89974e-06,9.15827e-10,-6.6401e-14,20959,-23.1458], Tmin=(985.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(ROOJ)"""),
)

species(
    label = '[CH2][CH]O[O](697)',
    structure = SMILES('[CH2][CH]O[O]'),
    E0 = (369.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.141847,'amu*angstrom^2'), symmetry=1, barrier=(3.26134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241494,'amu*angstrom^2'), symmetry=1, barrier=(5.55241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57149,0.0354185,-5.6472e-05,5.05292e-08,-1.74356e-11,44486.2,16.611], Tmin=(100,'K'), Tmax=(876.907,'K')), NASAPolynomial(coeffs=[5.01042,0.0153745,-6.92922e-06,1.2659e-09,-8.43558e-14,44401.4,7.12024], Tmin=(876.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CJCOOH) + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH]C=C(8168)',
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
    label = '[CH]OO[CH][CH2](1003)',
    structure = SMILES('[CH]OO[CH][CH2]'),
    E0 = (676.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,500,795,815,3000,3100,440,815,1455,1000,322.921,322.933,2148.64],'cm^-1')),
        HinderedRotor(inertia=(0.136318,'amu*angstrom^2'), symmetry=1, barrier=(10.0871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00987465,'amu*angstrom^2'), symmetry=1, barrier=(32.3493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.43713,'amu*angstrom^2'), symmetry=1, barrier=(32.3494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00987456,'amu*angstrom^2'), symmetry=1, barrier=(32.3497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66653,0.0559057,-8.48629e-05,7.19246e-08,-2.47308e-11,81413.7,21.3333], Tmin=(100,'K'), Tmax=(757.305,'K')), NASAPolynomial(coeffs=[7.71562,0.0212676,-1.09318e-05,2.15597e-09,-1.52002e-13,80574.6,-5.65971], Tmin=(757.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CJCOOH) + radical(CCsJOOC) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2][C]OO[CH]C=C(14698)',
    structure = SMILES('[CH2][C]OO[CH]C=C'),
    E0 = (679.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,350,500,795,815,3010,987.5,1337.5,450,1655,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.446967,0.0698142,-6.44129e-05,2.8781e-08,-5.03914e-12,81813.5,29.8549], Tmin=(100,'K'), Tmax=(1382.63,'K')), NASAPolynomial(coeffs=[18.48,0.0176445,-7.81498e-06,1.49127e-09,-1.04785e-13,76826.9,-62.9865], Tmin=(1382.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(CJCOOH) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2][CH]OO[C]=C[CH2](14699)',
    structure = SMILES('[CH2][CH]OO[C]=C[CH2]'),
    E0 = (675.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,350,500,795,815,1685,370,3010,987.5,1337.5,450,1655,1962.47],'cm^-1')),
        HinderedRotor(inertia=(0.640502,'amu*angstrom^2'), symmetry=1, barrier=(14.7264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.482049,'amu*angstrom^2'), symmetry=1, barrier=(11.0833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0158413,'amu*angstrom^2'), symmetry=1, barrier=(43.2928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286676,'amu*angstrom^2'), symmetry=1, barrier=(43.2926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88299,'amu*angstrom^2'), symmetry=1, barrier=(43.2936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05516,0.0698009,-8.29865e-05,5.71793e-08,-1.65177e-11,81387.4,31.4393], Tmin=(100,'K'), Tmax=(830.47,'K')), NASAPolynomial(coeffs=[8.8109,0.0324453,-1.55154e-05,3.01683e-09,-2.13101e-13,80099.2,-4.53702], Tmin=(830.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(675.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CJO) + radical(CJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = '[CH][CH]OO[CH]C=C(14700)',
    structure = SMILES('[CH][CH]OO[CH]C=C'),
    E0 = (638.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,350,500,795,815,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.836865,0.0642536,-5.36746e-05,2.22062e-08,-3.68184e-12,76950.4,29.0239], Tmin=(100,'K'), Tmax=(1428.14,'K')), NASAPolynomial(coeffs=[15.5894,0.0229343,-1.02765e-05,1.9477e-09,-1.35556e-13,72736.7,-47.406], Tmin=(1428.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(CCsJOOC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C1OOC1[CH2](14701)',
    structure = SMILES('[CH2][CH]C1OOC1[CH2]'),
    E0 = (513.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0256424,0.0669735,-5.93377e-05,2.75874e-08,-4.854e-12,61890.1,32.4407], Tmin=(100,'K'), Tmax=(1646.82,'K')), NASAPolynomial(coeffs=[15.2473,0.0175882,-3.04894e-06,2.23499e-10,-5.1193e-15,58559.9,-43.4775], Tmin=(1646.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxetane) + radical(CJCOOH) + radical(CCJCOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1C[CH]OO1(14702)',
    structure = SMILES('[CH2][CH]C1C[CH]OO1'),
    E0 = (407.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41618,0.0458509,4.85503e-06,-4.57385e-08,2.40879e-11,49111.7,26.0667], Tmin=(100,'K'), Tmax=(886.293,'K')), NASAPolynomial(coeffs=[13.826,0.0206383,-4.59335e-06,5.72538e-10,-3.36291e-14,45702.5,-39.1291], Tmin=(886.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CCJCOOH) + radical(RCCJ) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2][CH]OOC=[C]C(14703)',
    structure = SMILES('[CH2][CH]OOC=[C]C'),
    E0 = (522.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.812667,0.0751509,-8.95893e-05,6.13931e-08,-1.75257e-11,62946.2,29.3737], Tmin=(100,'K'), Tmax=(842.326,'K')), NASAPolynomial(coeffs=[9.52633,0.0337719,-1.5902e-05,3.0726e-09,-2.16296e-13,61478.3,-11.1696], Tmin=(842.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOOC) + radical(CJCOOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]OO[C]=CC(14704)',
    structure = SMILES('[CH2][CH]OO[C]=CC'),
    E0 = (524.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59722,0.0631037,-3.98961e-05,-4.45853e-08,6.51444e-11,63140.7,29.3599], Tmin=(100,'K'), Tmax=(482.329,'K')), NASAPolynomial(coeffs=[6.81691,0.0379756,-1.82236e-05,3.51611e-09,-2.45742e-13,62426,5.79405], Tmin=(482.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOOC) + radical(CJCOOH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=[C]OOC[CH2](14705)',
    structure = SMILES('[CH2]C=[C]OOC[CH2]'),
    E0 = (489.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,500,795,815,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.166807,'amu*angstrom^2'), symmetry=1, barrier=(3.83521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58072,'amu*angstrom^2'), symmetry=1, barrier=(36.3438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.01379,'amu*angstrom^2'), symmetry=1, barrier=(3.81763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.431985,'amu*angstrom^2'), symmetry=1, barrier=(9.9322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183363,'amu*angstrom^2'), symmetry=1, barrier=(51.4964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12306,0.0671187,-6.61487e-05,3.66719e-08,-8.60406e-12,58962.6,30.897], Tmin=(100,'K'), Tmax=(1003.19,'K')), NASAPolynomial(coeffs=[9.62179,0.033231,-1.54777e-05,2.99777e-09,-2.12098e-13,57257.5,-10.1313], Tmin=(1003.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CJCOOH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=[C]OO[CH]C(13631)',
    structure = SMILES('[CH2]C=[C]OO[CH]C'),
    E0 = (461.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,500,795,815,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1685,370,3010,987.5,1337.5,450,1655,860.883],'cm^-1')),
        HinderedRotor(inertia=(0.387277,'amu*angstrom^2'), symmetry=1, barrier=(8.90425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92071,'amu*angstrom^2'), symmetry=1, barrier=(44.1608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60586,'amu*angstrom^2'), symmetry=1, barrier=(13.9299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00806736,'amu*angstrom^2'), symmetry=1, barrier=(44.1597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92079,'amu*angstrom^2'), symmetry=1, barrier=(44.1628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24188,0.065625,-6.4608e-05,3.76112e-08,-9.53324e-12,55647.2,28.7759], Tmin=(100,'K'), Tmax=(921.102,'K')), NASAPolynomial(coeffs=[7.93227,0.036571,-1.72938e-05,3.36642e-09,-2.38675e-13,54414.7,-2.95148], Tmin=(921.102,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CCsJOOC) + radical(C=CJO)"""),
)

species(
    label = 'C1=COC1(5259)',
    structure = SMILES('C1=COC1'),
    E0 = (4.81952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10338,-0.00295546,9.94696e-05,-1.38902e-07,5.64804e-11,633.036,9.26427], Tmin=(100,'K'), Tmax=(918.619,'K')), NASAPolynomial(coeffs=[16.7387,-0.00374949,5.11332e-06,-1.00759e-09,6.07118e-14,-4343.72,-68.8138], Tmin=(918.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.81952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]C1CC=COO1(14706)',
    structure = SMILES('[CH2]C1CC=COO1'),
    E0 = (139.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53865,0.0368851,3.59749e-05,-7.97795e-08,3.56185e-11,16913.6,21.0771], Tmin=(100,'K'), Tmax=(925.777,'K')), NASAPolynomial(coeffs=[16.371,0.0171151,-3.79589e-06,5.66923e-10,-4.16459e-14,12268.3,-59.5932], Tmin=(925.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][CH]OC([CH2])C=O(12656)',
    structure = SMILES('[CH2][CH]OC([CH2])C=O'),
    E0 = (212.549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3599.51,'J/mol'), sigma=(6.21639,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.23 K, Pc=34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.344482,0.0813666,-9.69486e-05,5.9741e-08,-1.45887e-11,25694.6,30.332], Tmin=(100,'K'), Tmax=(1000.09,'K')), NASAPolynomial(coeffs=[14.8053,0.0235287,-1.01999e-05,1.91386e-09,-1.33231e-13,22802.2,-39.4344], Tmin=(1000.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CJC(C)OC) + radical(CCsJOCs)"""),
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
    label = '[CH]=COO[CH][CH2](4629)',
    structure = SMILES('[CH]=COO[CH][CH2]'),
    E0 = (567.727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,350,500,795,815,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.287637,'amu*angstrom^2'), symmetry=1, barrier=(6.61334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0913296,'amu*angstrom^2'), symmetry=1, barrier=(22.0635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30631,'amu*angstrom^2'), symmetry=1, barrier=(30.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.15366,'amu*angstrom^2'), symmetry=1, barrier=(72.509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3381,0.0595805,-6.52458e-05,3.73576e-08,-8.61089e-12,68376.9,25.6014], Tmin=(100,'K'), Tmax=(1047.88,'K')), NASAPolynomial(coeffs=[11.5683,0.0205299,-9.3472e-06,1.79504e-09,-1.26614e-13,66232.8,-24.2321], Tmin=(1047.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCsJOOC) + radical(CJCOOH)"""),
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
    E0 = (404.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (524.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (567.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (609.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (677.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (666.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (607.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (571.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (616.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (569.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (451.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (854.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (863.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (630.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (471.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (431.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (429.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (412.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (404.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (755.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (751.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (999.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (890.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (887.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (850.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (526.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (515.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (404.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (612.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (659.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (590.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (604.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (487.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (412.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (718.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (983.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['C=C[O](594)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2]C1[CH]OOC1[CH2](14687)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.40382e+09,'s^-1'), n=0.352, Ea=(120.148,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2]C1[CH]OO[CH]C1(14688)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.33723e+10,'s^-1'), n=0.316667, Ea=(162.479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C=[C][CH]OOC=C(14336)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1092.27,'m^3/(mol*s)'), n=1.64867, Ea=(13.1815,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]OOC[C]=C(14689)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.09427e+10,'s^-1'), n=1.04582, Ea=(152.506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=CCOO[CH][CH2](14690)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.65416e+07,'s^-1'), n=1.654, Ea=(131.736,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]COO[CH][C]=C(14691)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.62365e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C][CH]OO[CH]C(13626)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6Hall;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C[CH]OOC[CH2](14692)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.73886e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C[CH]OO[CH]C(13627)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH][O](719)', '[CH2]C=C[O](5266)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(8)', '[CH2][CH]OO[CH][C]=C(14693)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(8)', '[CH]=C[CH]OO[CH][CH2](14694)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2][CH]OOC1[CH]C1(12044)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_csHO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2]C1C[CH][CH]OO1(14695)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.52313e+07,'s^-1'), n=0.767814, Ea=(66.7931,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH]1[CH]OO[CH]CC1(14696)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2]COOC=C=C(14697)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2]C1OOC1C=C(12657)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['C=C[CH]OOC=C(12635)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH][CH2](721)', 'C=C[CH]O[O](6572)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]O[O](697)', '[CH]C=C(8168)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(64)', '[CH]OO[CH][CH2](1003)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2][C]OO[CH]C=C(14698)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH2][CH]OO[C]=C[CH2](14699)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH][CH]OO[CH]C=C(14700)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2][CH]C1OOC1[CH2](14701)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2][CH]C1C[CH]OO1(14702)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.13873e+09,'s^-1'), n=0.337103, Ea=(111.427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[O](594)', '[CH2]C=C[O](5266)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.3e+11,'cm^3/(mol*s)'), n=0, Ea=(339.456,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_R;O_rad/OneDe] for rate rule [Od_R;O_rad/OneDe]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 336.5 to 339.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2][CH]OOC=[C]C(14703)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2][CH]OO[C]=CC(14704)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C=[C]OOC[CH2](14705)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.572e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C=[C]OO[CH]C(13631)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2][CH][O](719)', 'C1=COC1(5259)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO_intra] for rate rule [R3OO_SD;C_pri_rad_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2]C1CC=COO1(14706)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using an average for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]OO[CH]C=C(12655)'],
    products = ['[CH2][CH]OC([CH2])C=O(12656)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(T)(28)', '[CH]=COO[CH][CH2](4629)'],
    products = ['[CH2][CH]OO[CH]C=C(12655)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '3408',
    isomers = [
        '[CH2][CH]OO[CH]C=C(12655)',
    ],
    reactants = [
        ('C=C[O](594)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3408',
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

