species(
    label = '[CH2]C([O])C(O)[C]=C(14793)',
    structure = SMILES('[CH2]C([O])C(O)[C]=C'),
    E0 = (294.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,427.651,427.885,428.207],'cm^-1')),
        HinderedRotor(inertia=(0.000922752,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916945,'amu*angstrom^2'), symmetry=1, barrier=(11.9114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916249,'amu*angstrom^2'), symmetry=1, barrier=(11.9094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0917068,'amu*angstrom^2'), symmetry=1, barrier=(11.9019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671625,0.071628,-7.41515e-05,3.99426e-08,-8.59586e-12,35506.8,30.812], Tmin=(100,'K'), Tmax=(1125.09,'K')), NASAPolynomial(coeffs=[14.2408,0.0233856,-9.83334e-06,1.83106e-09,-1.27278e-13,32453.5,-36.2509], Tmin=(1125.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S)"""),
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
    label = '[CH2]C([O])C(O)=C=C(28666)',
    structure = SMILES('[CH2]C([O])C(O)=C=C'),
    E0 = (178.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.813292,'amu*angstrom^2'), symmetry=1, barrier=(18.6992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813825,'amu*angstrom^2'), symmetry=1, barrier=(18.7114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812402,'amu*angstrom^2'), symmetry=1, barrier=(18.6787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.217506,0.0766229,-8.06511e-05,4.03232e-08,-7.57678e-12,21651.2,29.8498], Tmin=(100,'K'), Tmax=(1452.66,'K')), NASAPolynomial(coeffs=[21.9152,0.00770143,-1.24617e-06,1.01754e-10,-4.15449e-15,16062.7,-82.2951], Tmin=(1452.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'C=[C]C(O)C(=C)[O](14497)',
    structure = SMILES('C=[C]C(O)C(=C)[O]'),
    E0 = (72.0046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,418.502,418.788,419.084],'cm^-1')),
        HinderedRotor(inertia=(0.0748184,'amu*angstrom^2'), symmetry=1, barrier=(9.33062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0751913,'amu*angstrom^2'), symmetry=1, barrier=(9.33013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749943,'amu*angstrom^2'), symmetry=1, barrier=(9.32591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4228.96,'J/mol'), sigma=(6.78198,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=660.55 K, Pc=30.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.910647,0.0639176,-6.36516e-05,3.26891e-08,-6.65416e-12,8774.8,29.0222], Tmin=(100,'K'), Tmax=(1193.58,'K')), NASAPolynomial(coeffs=[14.2257,0.0192951,-7.57309e-06,1.36659e-09,-9.3486e-14,5596.3,-37.5713], Tmin=(1193.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.0046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#CC(O)C([CH2])[O](23351)',
    structure = SMILES('C#CC(O)C([CH2])[O]'),
    E0 = (222.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,750,770,3400,2100,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,328.355,328.731],'cm^-1')),
        HinderedRotor(inertia=(0.00156458,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228029,'amu*angstrom^2'), symmetry=1, barrier=(17.5131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45846,'amu*angstrom^2'), symmetry=1, barrier=(35.0671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226545,'amu*angstrom^2'), symmetry=1, barrier=(17.5189,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491027,0.0721781,-7.65305e-05,3.85975e-08,-6.71762e-12,26869.9,28.6003], Tmin=(100,'K'), Tmax=(946.109,'K')), NASAPolynomial(coeffs=[16.4585,0.0159039,-5.12157e-06,8.29856e-10,-5.41111e-14,23345.7,-50.2061], Tmin=(946.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = 'C=[C]C(O)C=O(27744)',
    structure = SMILES('C=[C]C(O)C=O'),
    E0 = (-32.4957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,616.43],'cm^-1')),
        HinderedRotor(inertia=(0.478592,'amu*angstrom^2'), symmetry=1, barrier=(11.0038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0407669,'amu*angstrom^2'), symmetry=1, barrier=(10.9798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409213,'amu*angstrom^2'), symmetry=1, barrier=(10.9972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78798,0.043671,-3.37416e-05,1.31659e-08,-2.06728e-12,-3824.66,24.676], Tmin=(100,'K'), Tmax=(1501.82,'K')), NASAPolynomial(coeffs=[11.867,0.0168258,-6.92882e-06,1.26348e-09,-8.59241e-14,-6852.02,-28.0484], Tmin=(1501.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.4957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = 'C=[C]C(O)C=C(27708)',
    structure = SMILES('C=[C]C(O)C=C'),
    E0 = (148.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,401.788,402.949],'cm^-1')),
        HinderedRotor(inertia=(0.102172,'amu*angstrom^2'), symmetry=1, barrier=(11.7321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102173,'amu*angstrom^2'), symmetry=1, barrier=(11.7215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101377,'amu*angstrom^2'), symmetry=1, barrier=(11.7041,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43816,0.050265,-3.66185e-05,1.35366e-08,-2.02373e-12,17988.7,25.5802], Tmin=(100,'K'), Tmax=(1564.23,'K')), NASAPolynomial(coeffs=[12.9464,0.0208365,-8.39847e-06,1.50935e-09,-1.0151e-13,14388.3,-35.0892], Tmin=(1564.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = '[CH2]C([O])C=C=C(17827)',
    structure = SMILES('[CH2]C([O])C=C=C'),
    E0 = (393.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,288.104,288.969],'cm^-1')),
        HinderedRotor(inertia=(0.297564,'amu*angstrom^2'), symmetry=1, barrier=(17.7701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300663,'amu*angstrom^2'), symmetry=1, barrier=(17.7804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29603,0.05195,-3.77204e-05,9.62971e-09,5.10053e-13,47394.8,24.3535], Tmin=(100,'K'), Tmax=(1065.15,'K')), NASAPolynomial(coeffs=[13.692,0.0177628,-6.9876e-06,1.29203e-09,-9.07815e-14,44052.8,-39.5248], Tmin=(1065.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=C(O)C([CH2])[O](14791)',
    structure = SMILES('[CH2]C=C(O)C([CH2])[O]'),
    E0 = (153.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303283,0.0680129,-4.10274e-05,-9.20461e-09,1.22509e-11,18611.4,29.1991], Tmin=(100,'K'), Tmax=(943.758,'K')), NASAPolynomial(coeffs=[21.4335,0.0111708,-2.68091e-06,4.38686e-10,-3.35936e-14,13166.1,-79.2376], Tmin=(943.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](O)C(O)[C]=C(28370)',
    structure = SMILES('[CH2][C](O)C(O)[C]=C'),
    E0 = (240.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,267.918,268.758],'cm^-1')),
        HinderedRotor(inertia=(0.00234911,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196221,'amu*angstrom^2'), symmetry=1, barrier=(10.0251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196738,'amu*angstrom^2'), symmetry=1, barrier=(10.0277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196773,'amu*angstrom^2'), symmetry=1, barrier=(10.026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00235089,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.356008,0.0842841,-0.00011667,8.65428e-08,-2.55627e-11,29050.4,31.5655], Tmin=(100,'K'), Tmax=(830.038,'K')), NASAPolynomial(coeffs=[12.4207,0.026145,-1.16059e-05,2.15984e-09,-1.47827e-13,27047.6,-24.3923], Tmin=(830.038,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)[C](C)[O](28173)',
    structure = SMILES('C=[C]C(O)[C](C)[O]'),
    E0 = (259.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,206.25,206.282,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00396257,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569591,'amu*angstrom^2'), symmetry=1, barrier=(17.1921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569342,'amu*angstrom^2'), symmetry=1, barrier=(17.1923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569547,'amu*angstrom^2'), symmetry=1, barrier=(17.1924,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.804348,0.0734637,-8.49639e-05,5.39971e-08,-1.39909e-11,31293.1,29.1353], Tmin=(100,'K'), Tmax=(932.629,'K')), NASAPolynomial(coeffs=[11.2604,0.0286184,-1.28373e-05,2.4396e-09,-1.70592e-13,29342.8,-20.5802], Tmin=(932.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = '[CH]=CC(O)C([CH2])[O](14794)',
    structure = SMILES('[CH]=CC(O)C([CH2])[O]'),
    E0 = (303.467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,388.052,388.054],'cm^-1')),
        HinderedRotor(inertia=(0.126583,'amu*angstrom^2'), symmetry=1, barrier=(13.5266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111948,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126585,'amu*angstrom^2'), symmetry=1, barrier=(13.5266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126584,'amu*angstrom^2'), symmetry=1, barrier=(13.5266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.451921,0.0726671,-7.34203e-05,3.75996e-08,-7.57933e-12,36630.9,31.5288], Tmin=(100,'K'), Tmax=(1209.09,'K')), NASAPolynomial(coeffs=[16.4436,0.0197629,-7.78797e-06,1.41172e-09,-9.69616e-14,32763.8,-48.6582], Tmin=(1209.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])C([O])C=C(12643)',
    structure = SMILES('[CH2]C([O])C([O])C=C'),
    E0 = (286.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,443.436,443.438,443.441,443.441],'cm^-1')),
        HinderedRotor(inertia=(0.000857287,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112652,'amu*angstrom^2'), symmetry=1, barrier=(15.7197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11266,'amu*angstrom^2'), symmetry=1, barrier=(15.7197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4232,'J/mol'), sigma=(7.00504,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=661.03 K, Pc=27.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582469,0.0663766,-5.50375e-05,2.04652e-08,-2.25608e-12,34616.4,30.4667], Tmin=(100,'K'), Tmax=(1099.34,'K')), NASAPolynomial(coeffs=[16.3838,0.0203947,-8.0051e-06,1.46921e-09,-1.02434e-13,30446.5,-50.426], Tmin=(1099.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][C]([O])C(O)C=C(14792)',
    structure = SMILES('[CH2][C]([O])C(O)C=C'),
    E0 = (232.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.617175,0.0749332,-8.3347e-05,4.87445e-08,-1.14087e-11,28144.7,29.9614], Tmin=(100,'K'), Tmax=(1037.12,'K')), NASAPolynomial(coeffs=[13.6556,0.0246458,-1.06152e-05,1.99178e-09,-1.38798e-13,25440.2,-33.4166], Tmin=(1037.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]=C(O)C([CH2])O(28371)',
    structure = SMILES('[CH2][C]=C(O)C([CH2])O'),
    E0 = (161.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3580,3650,1210,1345,900,1100,350,440,435,1725,1380,1390,370,380,2900,435,327.99],'cm^-1')),
        HinderedRotor(inertia=(0.219008,'amu*angstrom^2'), symmetry=1, barrier=(16.7177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00156724,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219055,'amu*angstrom^2'), symmetry=1, barrier=(16.7162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219067,'amu*angstrom^2'), symmetry=1, barrier=(16.7172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219021,'amu*angstrom^2'), symmetry=1, barrier=(16.7164,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.30247,0.0813913,-8.82175e-05,4.61482e-08,-9.15159e-12,19532.1,32.0413], Tmin=(100,'K'), Tmax=(1341.23,'K')), NASAPolynomial(coeffs=[21.696,0.0101404,-2.22006e-06,2.65152e-10,-1.4373e-14,14138.8,-78.6547], Tmin=(1341.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJCO) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=C(O)C(C)[O](28174)',
    structure = SMILES('[CH2][C]=C(O)C(C)[O]'),
    E0 = (179.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,422,422.15],'cm^-1')),
        HinderedRotor(inertia=(0.128499,'amu*angstrom^2'), symmetry=1, barrier=(16.2534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128506,'amu*angstrom^2'), symmetry=1, barrier=(16.2526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128573,'amu*angstrom^2'), symmetry=1, barrier=(16.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128435,'amu*angstrom^2'), symmetry=1, barrier=(16.2514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.418081,0.0674797,-4.62809e-05,1.21988e-09,7.28759e-12,21762.8,28.6268], Tmin=(100,'K'), Tmax=(959.405,'K')), NASAPolynomial(coeffs=[19.4749,0.0143755,-4.44916e-06,7.77322e-10,-5.62065e-14,16893.6,-68.8411], Tmin=(959.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)C([O])[C]=C(13927)',
    structure = SMILES('[CH2]C(O)C([O])[C]=C'),
    E0 = (294.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,427.651,427.885,428.207],'cm^-1')),
        HinderedRotor(inertia=(0.000922752,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916945,'amu*angstrom^2'), symmetry=1, barrier=(11.9114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916249,'amu*angstrom^2'), symmetry=1, barrier=(11.9094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0917068,'amu*angstrom^2'), symmetry=1, barrier=(11.9019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671625,0.071628,-7.41515e-05,3.99426e-08,-8.59586e-12,35506.8,30.812], Tmin=(100,'K'), Tmax=(1125.09,'K')), NASAPolynomial(coeffs=[14.2408,0.0233856,-9.83334e-06,1.83106e-09,-1.27278e-13,32453.5,-36.2509], Tmin=(1125.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = 'C=[C]C([O])C(C)[O](13760)',
    structure = SMILES('C=[C]C([O])C(C)[O]'),
    E0 = (312.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2950,3100,1380,975,1025,1650,409.73,410.078,410.529,410.772],'cm^-1')),
        HinderedRotor(inertia=(0.00100242,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114223,'amu*angstrom^2'), symmetry=1, barrier=(13.5956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113882,'amu*angstrom^2'), symmetry=1, barrier=(13.5949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803562,0.0645657,-5.5657e-05,2.45519e-08,-4.34299e-12,37763.1,29.5141], Tmin=(100,'K'), Tmax=(1350.85,'K')), NASAPolynomial(coeffs=[14.8711,0.0229105,-9.40281e-06,1.72473e-09,-1.18418e-13,33962.4,-42.5842], Tmin=(1350.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(O)C([CH2])O(28372)',
    structure = SMILES('[CH]=[C]C(O)C([CH2])O'),
    E0 = (310.948,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,3580,3650,1210,1345,900,1100,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,356.389],'cm^-1')),
        HinderedRotor(inertia=(0.00132719,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124248,'amu*angstrom^2'), symmetry=1, barrier=(11.1991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124249,'amu*angstrom^2'), symmetry=1, barrier=(11.1991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124236,'amu*angstrom^2'), symmetry=1, barrier=(11.1989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124241,'amu*angstrom^2'), symmetry=1, barrier=(11.1991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.477751,0.0785847,-9.45383e-05,5.93166e-08,-1.47588e-11,37524.3,32.1072], Tmin=(100,'K'), Tmax=(982.228,'K')), NASAPolynomial(coeffs=[14.0967,0.0231225,-9.83876e-06,1.82795e-09,-1.26435e-13,34849,-33.3523], Tmin=(982.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = '[CH]=[C]C(O)C(C)[O](28175)',
    structure = SMILES('[CH]=[C]C(O)C(C)[O]'),
    E0 = (329.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,336.784,340.604],'cm^-1')),
        HinderedRotor(inertia=(0.00146322,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14965,'amu*angstrom^2'), symmetry=1, barrier=(12.1507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149663,'amu*angstrom^2'), symmetry=1, barrier=(12.1101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151104,'amu*angstrom^2'), symmetry=1, barrier=(12.1176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769319,0.0697327,-7.02388e-05,3.70571e-08,-7.85705e-12,39773.6,30.2308], Tmin=(100,'K'), Tmax=(1137.17,'K')), NASAPolynomial(coeffs=[13.6837,0.0243061,-1.03182e-05,1.92863e-09,-1.34254e-13,36836.4,-33.7339], Tmin=(1137.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])[CH][C]=C(17829)',
    structure = SMILES('[CH2]C([O])[CH][C]=C'),
    E0 = (583.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,427.125,427.437,428.328],'cm^-1')),
        HinderedRotor(inertia=(0.190337,'amu*angstrom^2'), symmetry=1, barrier=(24.6372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190245,'amu*angstrom^2'), symmetry=1, barrier=(24.6462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919187,'amu*angstrom^2'), symmetry=1, barrier=(119.625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.28011,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11284,0.0614599,-6.0597e-05,3.09619e-08,-6.34264e-12,70336.9,23.1595], Tmin=(100,'K'), Tmax=(1176.36,'K')), NASAPolynomial(coeffs=[13.0376,0.0209119,-8.89332e-06,1.66036e-09,-1.15459e-13,67531.4,-36.3072], Tmin=(1176.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(583.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C([O])[C]=C(14512)',
    structure = SMILES('[CH2]C([O])C([O])[C]=C'),
    E0 = (524.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,449.847,449.901,449.908,449.972],'cm^-1')),
        HinderedRotor(inertia=(0.000832594,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970294,'amu*angstrom^2'), symmetry=1, barrier=(13.9379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097018,'amu*angstrom^2'), symmetry=1, barrier=(13.9387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.810303,0.0679251,-6.92787e-05,3.62448e-08,-7.55064e-12,63208.4,30.1217], Tmin=(100,'K'), Tmax=(1162.03,'K')), NASAPolynomial(coeffs=[14.3589,0.0212856,-9.0721e-06,1.70243e-09,-1.18899e-13,60059.7,-37.2767], Tmin=(1162.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJCO) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][C]=C(O)C([CH2])[O](28328)',
    structure = SMILES('[CH2][C]=C(O)C([CH2])[O]'),
    E0 = (391.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,324.377,324.38],'cm^-1')),
        HinderedRotor(inertia=(0.00160213,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262492,'amu*angstrom^2'), symmetry=1, barrier=(19.6002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262503,'amu*angstrom^2'), symmetry=1, barrier=(19.6002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262506,'amu*angstrom^2'), symmetry=1, barrier=(19.6002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.189098,0.0779505,-8.41137e-05,4.32781e-08,-8.40043e-12,47234.9,31.4446], Tmin=(100,'K'), Tmax=(1385.9,'K')), NASAPolynomial(coeffs=[21.705,0.00821695,-1.55783e-06,1.59523e-10,-7.88088e-15,41794.6,-79.0599], Tmin=(1385.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])C(O)[C]=C(28667)',
    structure = SMILES('[CH2][C]([O])C(O)[C]=C'),
    E0 = (470.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,268.743,268.75,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00233408,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311045,'amu*angstrom^2'), symmetry=1, barrier=(15.942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419988,'amu*angstrom^2'), symmetry=1, barrier=(21.5257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311053,'amu*angstrom^2'), symmetry=1, barrier=(15.9421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547209,0.0799303,-0.000109383,7.95329e-08,-2.30278e-11,56749.7,30.6891], Tmin=(100,'K'), Tmax=(846.005,'K')), NASAPolynomial(coeffs=[12.2707,0.0245022,-1.11091e-05,2.0938e-09,-1.44648e-13,54766.1,-23.9096], Tmin=(846.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = '[CH]=[C]C(O)C([CH2])[O](28668)',
    structure = SMILES('[CH]=[C]C(O)C([CH2])[O]'),
    E0 = (541.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,3615,1277.5,1000,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,387.724,387.725],'cm^-1')),
        HinderedRotor(inertia=(0.00112132,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115031,'amu*angstrom^2'), symmetry=1, barrier=(12.2714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115033,'amu*angstrom^2'), symmetry=1, barrier=(12.2714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115031,'amu*angstrom^2'), symmetry=1, barrier=(12.2715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.643143,0.0745689,-8.85906e-05,5.42678e-08,-1.31606e-11,65224.7,31.3211], Tmin=(100,'K'), Tmax=(1007.23,'K')), NASAPolynomial(coeffs=[14.0633,0.0212727,-9.21874e-06,1.73215e-09,-1.20744e-13,62521.3,-33.5199], Tmin=(1007.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(O)C(=C)O(28367)',
    structure = SMILES('C=[C]C(O)C(=C)O'),
    E0 = (-65.8002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3580,3650,1210,1345,900,1100,350,440,435,1725,1380,1390,370,380,2900,435,320.053,320.098],'cm^-1')),
        HinderedRotor(inertia=(0.1936,'amu*angstrom^2'), symmetry=1, barrier=(14.0735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193622,'amu*angstrom^2'), symmetry=1, barrier=(14.0737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19361,'amu*angstrom^2'), symmetry=1, barrier=(14.0736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193617,'amu*angstrom^2'), symmetry=1, barrier=(14.0736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482716,0.0686369,-5.88485e-05,1.93162e-08,-3.14913e-13,-7779.58,29.056], Tmin=(100,'K'), Tmax=(994.518,'K')), NASAPolynomial(coeffs=[17.7696,0.0162507,-5.69169e-06,1.01505e-09,-7.13259e-14,-12065.8,-58.5105], Tmin=(994.518,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.8002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)C(C)=O(27722)',
    structure = SMILES('C=[C]C(O)C(C)=O'),
    E0 = (-87.3665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,1685,370,1380,1390,370,380,2900,435,335.61,335.626],'cm^-1')),
        HinderedRotor(inertia=(0.00149658,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00149706,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116063,'amu*angstrom^2'), symmetry=1, barrier=(9.27398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116048,'amu*angstrom^2'), symmetry=1, barrier=(9.27393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26888,0.0574432,-4.53322e-05,1.85772e-08,-3.1259e-12,-10407.5,28.1243], Tmin=(100,'K'), Tmax=(1384.49,'K')), NASAPolynomial(coeffs=[12.2349,0.0257608,-1.10067e-05,2.04874e-09,-1.41343e-13,-13444,-28.3484], Tmin=(1384.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.3665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)C(O)=C=C(28368)',
    structure = SMILES('[CH2]C(O)C(O)=C=C'),
    E0 = (-51.7193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.709008,'amu*angstrom^2'), symmetry=1, barrier=(16.3015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710113,'amu*angstrom^2'), symmetry=1, barrier=(16.3269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710477,'amu*angstrom^2'), symmetry=1, barrier=(16.3353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710292,'amu*angstrom^2'), symmetry=1, barrier=(16.331,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.323652,0.0799886,-8.45338e-05,4.2957e-08,-8.24581e-12,-6051.83,30.4199], Tmin=(100,'K'), Tmax=(1406.58,'K')), NASAPolynomial(coeffs=[21.9491,0.009566,-1.87941e-06,2.01306e-10,-1.01868e-14,-11616.7,-82.1407], Tmin=(1406.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.7193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCO)"""),
)

species(
    label = 'C=C=C(O)C(C)[O](28168)',
    structure = SMILES('C=C=C(O)C(C)[O]'),
    E0 = (-32.9476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.797876,'amu*angstrom^2'), symmetry=1, barrier=(18.3447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799273,'amu*angstrom^2'), symmetry=1, barrier=(18.3769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.798423,'amu*angstrom^2'), symmetry=1, barrier=(18.3573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491931,0.0650016,-3.90332e-05,-6.33726e-09,9.95401e-12,-3825.36,26.6612], Tmin=(100,'K'), Tmax=(958.209,'K')), NASAPolynomial(coeffs=[19.5686,0.0140817,-4.27289e-06,7.52602e-10,-5.52773e-14,-8799.5,-71.4371], Tmin=(958.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.9476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC(O)C(=C)[O](14798)',
    structure = SMILES('C=CC(O)C(=C)[O]'),
    E0 = (-165.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865519,0.0602501,-4.21973e-05,7.89108e-09,2.32314e-12,-19825.1,28.7098], Tmin=(100,'K'), Tmax=(1013.66,'K')), NASAPolynomial(coeffs=[15.3531,0.0198897,-7.34596e-06,1.32888e-09,-9.30576e-14,-23625.8,-45.6403], Tmin=(1013.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'C=[C]C(O)C[CH][O](14866)',
    structure = SMILES('C=[C]C(O)C[CH][O]'),
    E0 = (270.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,346.865,351.741,354.895,2099.78],'cm^-1')),
        HinderedRotor(inertia=(0.001309,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111825,'amu*angstrom^2'), symmetry=1, barrier=(9.85425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11228,'amu*angstrom^2'), symmetry=1, barrier=(9.84379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115972,'amu*angstrom^2'), symmetry=1, barrier=(9.91461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729806,0.0789591,-0.000116324,1.00749e-07,-3.51421e-11,32693.8,30.9144], Tmin=(100,'K'), Tmax=(813.682,'K')), NASAPolynomial(coeffs=[7.33732,0.0357866,-1.70293e-05,3.24803e-09,-2.24309e-13,31972.4,2.57398], Tmin=(813.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)C1CO1(27719)',
    structure = SMILES('C=[C]C(O)C1CO1'),
    E0 = (41.8263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24128,0.0661336,-6.06004e-05,2.95219e-08,-5.46684e-12,5180.06,28.9545], Tmin=(100,'K'), Tmax=(1551.2,'K')), NASAPolynomial(coeffs=[14.696,0.0177907,-3.14927e-06,2.30616e-10,-4.69029e-15,2027.38,-42.8345], Tmin=(1551.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.8263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1OC(=C)C1O(28652)',
    structure = SMILES('[CH2]C1OC(=C)C1O'),
    E0 = (-53.9572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.70352,0.0542631,3.89987e-06,-6.36784e-08,3.48679e-11,-6353.67,20.976], Tmin=(100,'K'), Tmax=(897.309,'K')), NASAPolynomial(coeffs=[22.9437,0.00506061,2.66796e-06,-7.39214e-10,5.18129e-14,-12355.4,-95.1132], Tmin=(897.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.9572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=C1CC([O])C1O(28669)',
    structure = SMILES('C=C1CC([O])C1O'),
    E0 = (-14.5626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.8141,-0.0151151,0.000127666,-1.23242e-07,2.84196e-11,-2114.76,-17.9819], Tmin=(100,'K'), Tmax=(1719.23,'K')), NASAPolynomial(coeffs=[76.1232,0.0217841,-6.92948e-05,1.70255e-08,-1.26791e-12,-52480.6,-449.619], Tmin=(1719.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.5626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([CH][C]=C)OO(28670)',
    structure = SMILES('[CH2]C([CH][C]=C)OO'),
    E0 = (427.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,227.245,227.246],'cm^-1')),
        HinderedRotor(inertia=(0.919383,'amu*angstrom^2'), symmetry=1, barrier=(33.6845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919337,'amu*angstrom^2'), symmetry=1, barrier=(33.6846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919209,'amu*angstrom^2'), symmetry=1, barrier=(33.6844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00326463,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919254,'amu*angstrom^2'), symmetry=1, barrier=(33.6845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.549935,0.077591,-8.25208e-05,4.61399e-08,-1.0456e-11,51533.7,27.5022], Tmin=(100,'K'), Tmax=(1060.67,'K')), NASAPolynomial(coeffs=[13.557,0.0285394,-1.31528e-05,2.54026e-09,-1.79706e-13,48774.4,-36.0156], Tmin=(1060.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJCO) + radical(CJCOOH)"""),
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
    label = '[CH2]C([O])[CH]O(1035)',
    structure = SMILES('[CH2]C([O])[CH]O'),
    E0 = (169.426,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3615,1277.5,1000,1380,1390,370,380,2900,435,383.027,383.164],'cm^-1')),
        HinderedRotor(inertia=(0.0011488,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0953298,'amu*angstrom^2'), symmetry=1, barrier=(9.92751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0953361,'amu*angstrom^2'), symmetry=1, barrier=(9.92798,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53775,0.057242,-7.71793e-05,5.29793e-08,-1.34641e-11,20463.1,21.2621], Tmin=(100,'K'), Tmax=(722.247,'K')), NASAPolynomial(coeffs=[10.1,0.0166363,-6.99961e-06,1.26347e-09,-8.48058e-14,19048.6,-18.4901], Tmin=(722.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CJCO) + radical(CC(C)OJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C(O)[C]=C(27815)',
    structure = SMILES('[CH2][CH]C(O)[C]=C'),
    E0 = (429.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,417.3,417.304],'cm^-1')),
        HinderedRotor(inertia=(0.000968012,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000967935,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000967927,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0703548,'amu*angstrom^2'), symmetry=1, barrier=(8.69422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25612,0.0525465,-4.15588e-05,1.68717e-08,-2.74125e-12,51808.1,30.2053], Tmin=(100,'K'), Tmax=(1469.61,'K')), NASAPolynomial(coeffs=[13.6561,0.018796,-7.1102e-06,1.2445e-09,-8.28526e-14,48163.5,-34.3914], Tmin=(1469.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)[CH][O](28102)',
    structure = SMILES('C=[C]C(O)[CH][O]'),
    E0 = (294.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,276.259,279.06,2001.85],'cm^-1')),
        HinderedRotor(inertia=(0.00227246,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191772,'amu*angstrom^2'), symmetry=1, barrier=(10.9158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189866,'amu*angstrom^2'), symmetry=1, barrier=(10.8672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40118,0.0638678,-0.000101269,9.07849e-08,-3.19727e-11,35530,26.2595], Tmin=(100,'K'), Tmax=(835.883,'K')), NASAPolynomial(coeffs=[6.37686,0.0274382,-1.32505e-05,2.52384e-09,-1.73332e-13,35139,5.78378], Tmin=(835.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S) + radical(CCsJOH)"""),
)

species(
    label = '[CH]C([O])C(O)[C]=C(28671)',
    structure = SMILES('[CH]C([O])C(O)[C]=C'),
    E0 = (530.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773224,0.0694706,-7.38227e-05,4.03262e-08,-8.75647e-12,63962.6,30.8179], Tmin=(100,'K'), Tmax=(1118.1,'K')), NASAPolynomial(coeffs=[14.2736,0.0211728,-9.02787e-06,1.69216e-09,-1.18096e-13,60943.7,-35.8208], Tmin=(1118.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCJ2_triplet) + radical(CC(C)OJ)"""),
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
    E0 = (294.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (402.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (305.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (447.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (404.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (391.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (357.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (294.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (421.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (476.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (408.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (435.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (408.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (436.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (436.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (369.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (412.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (399.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (377.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (343.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (362.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (612.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (743.473,'kJ/mol'),
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
    E0 = (603.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (682.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (753.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (317.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (317.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (372.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (372.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (383.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (388.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (451.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (299.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (302.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (302.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (534.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (803.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (836.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (710.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (742.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=C[O](594)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH2]C([O])C(O)=C=C(28666)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C=[C]C(O)C(=C)[O](14497)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC(O)C([CH2])[O](23351)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(28)', 'C=[C]C(O)C=O(27744)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(T)(63)', 'C=[C]C(O)C=C(27708)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;O_atom_triplet] for rate rule [Cds-Cs\O2s/H_Cds-HH;O_atom_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH][O](719)', 'C=C=CO(12571)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00401797,'m^3/(mol*s)'), n=2.41733, Ea=(22.1495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[O](594)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(132.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 129.1 to 132.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(D)(132)', '[CH2]C([O])C=C=C(17827)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(986500,'cm^3/(mol*s)'), n=2.037, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;OJ_pri] for rate rule [Cds-CsH_Ca;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2]C=C(O)C([CH2])[O](14791)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2][C](O)C(O)[C]=C(28370)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=[C]C(O)[C](C)[O](28173)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC(O)C([CH2])[O](14794)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2]C([O])C([O])C=C(12643)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2][C]([O])C(O)C=C(14792)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2][C]=C(O)C([CH2])O(28371)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2][C]=C(O)C(C)[O](28174)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)C([O])[C]=C(13927)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=[C]C([O])C(C)[O](13760)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]C(O)C([CH2])O(28372)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]C(O)C(C)[O](28175)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['OH(D)(132)', '[CH2]C([O])[CH][C]=C(17829)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2]C([O])C([O])[C]=C(14512)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.00518e+06,'m^3/(mol*s)'), n=0.282325, Ea=(7.09479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH][O](719)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH2][C]=C(O)C([CH2])[O](28328)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH2][C]([O])C(O)[C]=C(28667)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', '[CH]=[C]C(O)C([CH2])[O](28668)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=[C]C(O)C(=C)O(28367)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=[C]C(O)C(C)=O(27722)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2]C(O)C(O)=C=C(28368)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=C=C(O)C(C)[O](28168)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=CC(O)C(=C)[O](14798)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2]C(=CO)C([CH2])[O](14760)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=[C]C(O)C[CH][O](14866)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=[C]C(O)C1CO1(27719)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['[CH2]C1OC(=C)C1O(28652)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([O])C(O)[C]=C(14793)'],
    products = ['C=C1CC([O])C1O(28669)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH][C]=C)OO(28670)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[C]=C(584)', '[CH2]C([O])[CH]O(1035)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['O(T)(63)', '[CH2][CH]C(O)[C]=C(27815)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(T)(28)', 'C=[C]C(O)[CH][O](28102)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(8)', '[CH]C([O])C(O)[C]=C(28671)'],
    products = ['[CH2]C([O])C(O)[C]=C(14793)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '5035',
    isomers = [
        '[CH2]C([O])C(O)[C]=C(14793)',
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
    label = '5035',
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

