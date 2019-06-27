species(
    label = 'C=CC([O])C([O])=O(12739)',
    structure = SMILES('C=CC([O])C([O])=O'),
    E0 = (-66.5857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180,180,827.09,827.182,827.199,827.204,827.214],'cm^-1')),
        HinderedRotor(inertia=(0.00486802,'amu*angstrom^2'), symmetry=1, barrier=(2.36312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.004866,'amu*angstrom^2'), symmetry=1, barrier=(2.36278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33652,0.0389531,-2.20131e-05,5.44484e-09,-5.19577e-13,-7951.45,24.4788], Tmin=(100,'K'), Tmax=(2372.1,'K')), NASAPolynomial(coeffs=[16.0297,0.0158628,-7.41191e-06,1.34126e-09,-8.70945e-14,-14447.8,-53.4107], Tmin=(2372.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.5857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=OCOJ)"""),
)

species(
    label = 'O=C=O(1731)',
    structure = SMILES('O=C=O'),
    E0 = (-403.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.923,1087.69,1087.69,2296.71],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2028.74,'J/mol'), sigma=(3.763,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(2.65,'angstroms^3'), rotrelaxcollnum=2.1, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27862,0.0027414,7.16119e-06,-1.08033e-08,4.14308e-12,-48470.3,5.97933], Tmin=(100,'K'), Tmax=(988.876,'K')), NASAPolynomial(coeffs=[4.54605,0.0029192,-1.15488e-06,2.27663e-10,-1.70918e-14,-48980.3,-1.43251], Tmin=(988.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.087,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cdd-OdOd)"""),
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
    label = '[CH2]C1OC1C([O])=O(14648)',
    structure = SMILES('[CH2]C1OC1C([O])=O'),
    E0 = (-48.4564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88733,0.0374424,2.58814e-06,-3.48029e-08,1.82262e-11,-5743.43,22.1466], Tmin=(100,'K'), Tmax=(899.304,'K')), NASAPolynomial(coeffs=[11.9705,0.0170902,-4.324e-06,6.10372e-10,-3.85801e-14,-8547.56,-30.9359], Tmin=(899.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.4564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Ethylene_oxide) + radical(CJCO) + radical(CCOJ)"""),
)

species(
    label = 'C=CC1OC1([O])[O](14609)',
    structure = SMILES('C=CC1OC1([O])[O]'),
    E0 = (71.0365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92885,0.0463583,-4.75867e-05,3.24532e-08,-9.1827e-12,8617.48,23.8635], Tmin=(100,'K'), Tmax=(990.003,'K')), NASAPolynomial(coeffs=[5.82969,0.0257127,-8.9043e-06,1.4206e-09,-8.76437e-14,8084.49,6.29234], Tmin=(990.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.0365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1OC(=O)C1[O](14528)',
    structure = SMILES('[CH2]C1OC(=O)C1[O]'),
    E0 = (-48.5264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08723,0.025477,4.48814e-05,-8.70959e-08,3.90802e-11,-5751.94,24.2318], Tmin=(100,'K'), Tmax=(901.68,'K')), NASAPolynomial(coeffs=[15.9361,0.00803258,7.19054e-07,-3.36039e-10,2.31897e-14,-10037.7,-51.0642], Tmin=(901.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.5264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(C=OCOJ) + radical(CJC(C)OC)"""),
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
    label = 'C=CC(=O)C([O])=O(14649)',
    structure = SMILES('C=CC(=O)C([O])=O'),
    E0 = (-180.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,215.75,215.75,215.751,215.751,215.751,215.752,1121.9],'cm^-1')),
        HinderedRotor(inertia=(0.000133934,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0148168,'amu*angstrom^2'), symmetry=1, barrier=(13.234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76092,0.0545085,-7.12343e-05,5.56127e-08,-1.85171e-11,-21577.8,22.348], Tmin=(100,'K'), Tmax=(719.839,'K')), NASAPolynomial(coeffs=[6.6772,0.0271904,-1.43103e-05,2.89487e-09,-2.087e-13,-22285.6,0.245767], Tmin=(719.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(C=OC=OOJ)"""),
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
    label = '[O]C(=O)C=O(3180)',
    structure = SMILES('[O]C(=O)C=O'),
    E0 = (-220.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,556.455,679.051,679.566,3886.04],'cm^-1')),
        HinderedRotor(inertia=(1.37248,'amu*angstrom^2'), symmetry=1, barrier=(31.5561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0275,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.94679,0.0265977,-2.97427e-05,1.9056e-08,-5.37282e-12,-26451.4,14.3847], Tmin=(100,'K'), Tmax=(825.858,'K')), NASAPolynomial(coeffs=[5.38714,0.0147781,-8.27524e-06,1.72678e-09,-1.27056e-13,-26854.5,3.07836], Tmin=(825.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(C=OC=OOJ)"""),
)

species(
    label = '[O][C]=O(2059)',
    structure = SMILES('[O][C]=O'),
    E0 = (33.3014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81048,-0.00025715,1.76446e-05,-2.38747e-08,9.15883e-12,4016.03,8.55818], Tmin=(100,'K'), Tmax=(975.962,'K')), NASAPolynomial(coeffs=[6.50409,-1.44217e-05,-6.90664e-08,7.0435e-11,-9.1126e-15,2952.93,-7.12421], Tmin=(975.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-OdOsH) + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = 'C=CC(O)=C([O])[O](14650)',
    structure = SMILES('C=CC(O)=C([O])[O]'),
    E0 = (-202.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.025235,0.0737993,-8.48222e-05,4.50812e-08,-8.88172e-12,-24254.8,24.1535], Tmin=(100,'K'), Tmax=(1410.34,'K')), NASAPolynomial(coeffs=[21.6136,0.00278604,9.59977e-07,-3.15229e-10,2.45812e-14,-29399.5,-84.282], Tmin=(1410.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-202.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C(O)C([O])=O(14651)',
    structure = SMILES('C=[C]C(O)C([O])=O'),
    E0 = (-72.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25899,0.0439993,-3.39358e-05,1.47635e-08,-2.9929e-12,-8659.07,24.4669], Tmin=(100,'K'), Tmax=(1053.32,'K')), NASAPolynomial(coeffs=[5.6704,0.0310445,-1.54874e-05,3.0872e-09,-2.21601e-13,-9377.73,7.83167], Tmin=(1053.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC([O])=C([O])O(14652)',
    structure = SMILES('C=CC([O])=C([O])O'),
    E0 = (-206.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.029785,0.0757232,-9.2517e-05,5.25469e-08,-1.10616e-11,-24699.3,24.0399], Tmin=(100,'K'), Tmax=(1317.76,'K')), NASAPolynomial(coeffs=[20.8508,0.00293375,1.253e-06,-4.13728e-10,3.34243e-14,-29354.2,-78.9959], Tmin=(1317.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C([O])C(=O)O(14653)',
    structure = SMILES('C=[C]C([O])C(=O)O'),
    E0 = (-54.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83381,0.0475899,-3.94072e-05,1.69208e-08,-3.00618e-12,-6470.85,26.592], Tmin=(100,'K'), Tmax=(1306.07,'K')), NASAPolynomial(coeffs=[10.1895,0.0219995,-1.0017e-05,1.91889e-09,-1.34619e-13,-8653.49,-15.9507], Tmin=(1306.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(O)C([O])=O(14654)',
    structure = SMILES('[CH]=CC(O)C([O])=O'),
    E0 = (-63.2226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,242.407,242.498,242.728,242.86,1434.26,1434.44],'cm^-1')),
        HinderedRotor(inertia=(0.00286408,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0541333,'amu*angstrom^2'), symmetry=1, barrier=(35.6162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0243993,'amu*angstrom^2'), symmetry=1, barrier=(35.6149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13101,0.0442222,-3.13218e-05,1.09203e-08,-1.59645e-12,-7539.89,24.8333], Tmin=(100,'K'), Tmax=(1497.6,'K')), NASAPolynomial(coeffs=[9.3347,0.0249815,-1.20502e-05,2.34141e-09,-1.64337e-13,-9697.53,-12.8296], Tmin=(1497.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.2226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCOJ)"""),
)

species(
    label = '[CH]=CC([O])C(=O)O(14655)',
    structure = SMILES('[CH]=CC([O])C(=O)O'),
    E0 = (-45.1947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54183,0.0493853,-4.09095e-05,1.69637e-08,-2.81854e-12,-5343.4,27.575], Tmin=(100,'K'), Tmax=(1428.49,'K')), NASAPolynomial(coeffs=[12.8702,0.017664,-7.60031e-06,1.41853e-09,-9.79901e-14,-8579.91,-31.1179], Tmin=(1428.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.1947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_P)"""),
)

species(
    label = '[O]C=C([O])[O](3183)',
    structure = SMILES('[O]C=C([O])[O]'),
    E0 = (-108.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,180,180,180],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0275,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23477,0.0289406,-7.48908e-06,-2.38397e-08,1.47722e-11,-12938.9,15.4323], Tmin=(100,'K'), Tmax=(918.148,'K')), NASAPolynomial(coeffs=[15.6525,-0.00349641,2.99678e-06,-5.89033e-10,3.7326e-14,-16499.4,-54.1269], Tmin=(918.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C=CC([O])=C([O])[O](14656)',
    structure = SMILES('C=CC([O])=C([O])[O]'),
    E0 = (-65.171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.946809,'amu*angstrom^2'), symmetry=1, barrier=(21.769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.864539,0.0636848,-7.10181e-05,3.47393e-08,-5.30268e-12,-7720.53,22.4597], Tmin=(100,'K'), Tmax=(937.092,'K')), NASAPolynomial(coeffs=[16.9789,0.00774945,-2.05005e-06,3.06181e-10,-2.01427e-14,-11304.8,-57.246], Tmin=(937.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C([O])C([O])=O(14657)',
    structure = SMILES('C=[C]C([O])C([O])=O'),
    E0 = (171.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,180,180,180,905.244,905.274,905.473,905.743,906.154],'cm^-1')),
        HinderedRotor(inertia=(0.00793206,'amu*angstrom^2'), symmetry=1, barrier=(4.61502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00794015,'amu*angstrom^2'), symmetry=1, barrier=(4.62076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11059,0.045872,-5.57759e-05,4.88832e-08,-1.92402e-11,20661.3,25.7678], Tmin=(100,'K'), Tmax=(702.55,'K')), NASAPolynomial(coeffs=[4.08421,0.0308496,-1.56193e-05,3.10788e-09,-2.21943e-13,20477.4,17.6079], Tmin=(702.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH]=CC([O])C([O])=O(14658)',
    structure = SMILES('[CH]=CC([O])C([O])=O'),
    E0 = (180.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180,839.684,840.792,842.383,842.693,843.33],'cm^-1')),
        HinderedRotor(inertia=(0.00384592,'amu*angstrom^2'), symmetry=1, barrier=(1.93168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0865483,'amu*angstrom^2'), symmetry=1, barrier=(1.98992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.36301,0.0410952,-3.32998e-05,1.56513e-08,-3.44057e-12,21765.1,24.8059], Tmin=(100,'K'), Tmax=(987.945,'K')), NASAPolynomial(coeffs=[5.44378,0.0286217,-1.43611e-05,2.87137e-09,-2.06589e-13,21156.4,9.98038], Tmin=(987.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_P) + radical(CCOJ)"""),
)

species(
    label = 'C=CC([O])[C]1OO1(14659)',
    structure = SMILES('C=CC([O])[C]1OO1'),
    E0 = (280.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46663,0.0403294,9.07637e-06,-5.32793e-08,2.70483e-11,33833.6,26.346], Tmin=(100,'K'), Tmax=(923.458,'K')), NASAPolynomial(coeffs=[18.6032,0.00595568,1.7467e-07,-1.18462e-10,4.46432e-15,28969.3,-64.1642], Tmin=(923.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(dioxirane) + radical(CC(C)OJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C(=O)C1[CH]CO1(14660)',
    structure = SMILES('[O]C(=O)C1[CH]CO1'),
    E0 = (-53.0324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64072,0.0162351,5.66534e-05,-8.55117e-08,3.47687e-11,-6316.72,22.711], Tmin=(100,'K'), Tmax=(918.227,'K')), NASAPolynomial(coeffs=[10.1331,0.0187127,-4.75931e-06,7.25515e-10,-5.01308e-14,-9173.06,-20.8575], Tmin=(918.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.0324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Oxetane) + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = '[CH2][CH]C1OOC1=O(14553)',
    structure = SMILES('[CH2][CH]C1OOC1=O'),
    E0 = (68.7704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1871,0.0133205,9.57915e-05,-1.40703e-07,5.58234e-11,8360.55,27.0483], Tmin=(100,'K'), Tmax=(958.901,'K')), NASAPolynomial(coeffs=[19.3393,0.00893038,-2.39807e-06,6.02609e-10,-5.96588e-14,1983.47,-71.0813], Tmin=(958.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.7704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclobutane) + radical(CCJCOOH) + radical(RCCJ)"""),
)

species(
    label = '[O]C1[CH]COC1=O(14569)',
    structure = SMILES('[O]C1[CH]COC1=O'),
    E0 = (-108.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91006,0.000842405,0.000111748,-1.52705e-07,6.12711e-11,-12970.8,24.7099], Tmin=(100,'K'), Tmax=(913.723,'K')), NASAPolynomial(coeffs=[14.7582,0.00866909,9.0353e-07,-3.31067e-10,1.80619e-14,-17627.9,-45.0176], Tmin=(913.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(butyrolactone) + radical(C=OCOJ) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(=O)C(=O)O(14661)',
    structure = SMILES('C=CC(=O)C(=O)O'),
    E0 = (-441.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04557,0.0484065,-4.02145e-05,1.807e-08,-3.58531e-12,-53049.7,20.2379], Tmin=(100,'K'), Tmax=(1119.36,'K')), NASAPolynomial(coeffs=[7.33847,0.0294923,-1.48683e-05,2.97417e-09,-2.1376e-13,-54234.6,-5.89403], Tmin=(1119.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-441.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]CC([O])=C([O])[O](12041)',
    structure = SMILES('[CH2]CC([O])=C([O])[O]'),
    E0 = (28.9667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,487.101,487.104,487.121,4000],'cm^-1')),
        HinderedRotor(inertia=(0.110686,'amu*angstrom^2'), symmetry=1, barrier=(18.6412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110687,'amu*angstrom^2'), symmetry=1, barrier=(18.6406,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.936855,0.0634578,-7.12659e-05,3.99307e-08,-8.69505e-12,3597.56,26.8843], Tmin=(100,'K'), Tmax=(1130.14,'K')), NASAPolynomial(coeffs=[15.0223,0.0136048,-5.09862e-06,8.99372e-10,-6.09907e-14,413.815,-42.7932], Tmin=(1130.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.9667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C([O])C([O])[O](14662)',
    structure = SMILES('[CH2]C=C([O])C([O])[O]'),
    E0 = (136.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,180,1366.23,1366.47,1366.91],'cm^-1')),
        HinderedRotor(inertia=(0.276647,'amu*angstrom^2'), symmetry=1, barrier=(6.36067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276663,'amu*angstrom^2'), symmetry=1, barrier=(6.36103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42165,0.0645207,-0.000103036,9.63513e-08,-3.54574e-11,16465.8,26.9544], Tmin=(100,'K'), Tmax=(823.608,'K')), NASAPolynomial(coeffs=[4.86625,0.0322849,-1.60852e-05,3.10953e-09,-2.15626e-13,16424.3,14.1975], Tmin=(823.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCOJ) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C([O])C([O])[O](14663)',
    structure = SMILES('C=[C]C([O])C([O])[O]'),
    E0 = (370.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,1685,370,280.798,281.132,281.234,1452.92,1453.1],'cm^-1')),
        HinderedRotor(inertia=(0.00472749,'amu*angstrom^2'), symmetry=1, barrier=(7.0834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126194,'amu*angstrom^2'), symmetry=1, barrier=(7.08213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5166,0.0634386,-0.000104437,1.02694e-07,-3.96189e-11,44608.7,29.4272], Tmin=(100,'K'), Tmax=(807.41,'K')), NASAPolynomial(coeffs=[3.46295,0.035488,-1.84977e-05,3.65097e-09,-2.56404e-13,44891.2,24.1491], Tmin=(807.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]C([O])[C]([O])O(14664)',
    structure = SMILES('C=[C]C([O])[C]([O])O'),
    E0 = (349.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,360,370,350,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,254.783,254.799,1581.6,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0122437,'amu*angstrom^2'), symmetry=1, barrier=(21.7335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235461,'amu*angstrom^2'), symmetry=1, barrier=(10.8461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471781,'amu*angstrom^2'), symmetry=1, barrier=(21.7334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2445,0.0662286,-9.98358e-05,8.57169e-08,-2.98961e-11,42161.1,30.9286], Tmin=(100,'K'), Tmax=(767.502,'K')), NASAPolynomial(coeffs=[7.86106,0.0269351,-1.36403e-05,2.68039e-09,-1.8865e-13,41287.1,1.68127], Tmin=(767.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(Cs_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([O])C([O])[O](14665)',
    structure = SMILES('[CH]=CC([O])C([O])[O]'),
    E0 = (379.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,180,180,1428.06,1428.61],'cm^-1')),
        HinderedRotor(inertia=(0.12639,'amu*angstrom^2'), symmetry=1, barrier=(2.90595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127341,'amu*angstrom^2'), symmetry=1, barrier=(2.92781,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49435,0.0621297,-9.54487e-05,8.96807e-08,-3.4121e-11,45724.3,29.4376], Tmin=(100,'K'), Tmax=(779.907,'K')), NASAPolynomial(coeffs=[4.80934,0.0333066,-1.72777e-05,3.42559e-09,-2.42113e-13,45566.8,16.5735], Tmin=(779.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(379.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCOJ) + radical(CC(C)OJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]=CC([O])[C]([O])O(14666)',
    structure = SMILES('[CH]=CC([O])[C]([O])O'),
    E0 = (359.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,360,370,350,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,181.859,188.153,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139822,'amu*angstrom^2'), symmetry=1, barrier=(21.5038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.940323,'amu*angstrom^2'), symmetry=1, barrier=(21.6199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3577,0.0631134,-8.34044e-05,6.1084e-08,-1.83614e-11,43271,30.4657], Tmin=(100,'K'), Tmax=(805.385,'K')), NASAPolynomial(coeffs=[9.00726,0.0251207,-1.26433e-05,2.50971e-09,-1.79025e-13,42038.8,-4.78328], Tmin=(805.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cs_P) + radical(CCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC1OOC1=O(14551)',
    structure = SMILES('C=CC1OOC1=O'),
    E0 = (-210.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37472,0.0064117,0.000117677,-1.6445e-07,6.44707e-11,-25201.2,22.0895], Tmin=(100,'K'), Tmax=(955.41,'K')), NASAPolynomial(coeffs=[19.8338,0.00786086,-1.63343e-06,4.68103e-10,-5.21195e-14,-31939.6,-79.149], Tmin=(955.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=COC([O])=C[O](14667)',
    structure = SMILES('C=COC([O])=C[O]'),
    E0 = (-113.19,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14879,'amu*angstrom^2'), symmetry=1, barrier=(26.4129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14893,'amu*angstrom^2'), symmetry=1, barrier=(26.4163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.891702,0.0574313,-3.8312e-05,-3.99672e-09,9.10706e-12,-13491.7,24.6626], Tmin=(100,'K'), Tmax=(943.545,'K')), NASAPolynomial(coeffs=[19.0234,0.00744888,-1.59141e-06,2.55579e-10,-2.06658e-14,-18110.1,-68.1006], Tmin=(943.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-113.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C]1OC[CH]C1[O](14668)',
    structure = SMILES('[O][C]1OC[CH]C1[O]'),
    E0 = (289.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30797,0.0281188,1.96064e-05,-4.7493e-08,2.17414e-11,34942.5,22.9623], Tmin=(100,'K'), Tmax=(898.565,'K')), NASAPolynomial(coeffs=[9.96604,0.0183562,-4.70735e-06,6.75992e-10,-4.30716e-14,32584.1,-18.6294], Tmin=(898.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Tetrahydrofuran) + radical(Cs_P) + radical(CCJCO) + radical(CCOJ) + radical(CC(C)OJ)"""),
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
    label = 'C=C[CH]C([O])=O(5461)',
    structure = SMILES('C=C[CH]C([O])=O'),
    E0 = (-16.2033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,415.658,417.503,419.703,423.777,423.863,423.99,424.785],'cm^-1')),
        HinderedRotor(inertia=(0.000944513,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00094908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3065,0.0312387,-5.37794e-07,-1.43636e-08,5.12491e-12,-1883.14,17.1412], Tmin=(100,'K'), Tmax=(1302.84,'K')), NASAPolynomial(coeffs=[10.7396,0.0235631,-1.26731e-05,2.5778e-09,-1.85256e-13,-5626.51,-31.7077], Tmin=(1302.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.2033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(CCOJ)"""),
)

species(
    label = 'C=CC([O])[C]=O(14526)',
    structure = SMILES('C=CC([O])[C]=O'),
    E0 = (133.356,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,351.282,4000],'cm^-1')),
        HinderedRotor(inertia=(0.170973,'amu*angstrom^2'), symmetry=1, barrier=(14.9754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170895,'amu*angstrom^2'), symmetry=1, barrier=(14.9739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96577,0.0394789,-2.95633e-05,9.34197e-09,-5.31118e-13,16116.8,25.1576], Tmin=(100,'K'), Tmax=(1090.53,'K')), NASAPolynomial(coeffs=[10.8044,0.0145668,-5.62355e-06,1.01981e-09,-7.04881e-14,13742.6,-20.2963], Tmin=(1090.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
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
    E0 = (-66.5857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (41.5703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (71.4332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (55.4573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (53.3348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (81.7957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-28.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-66.5857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (58.8744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (109.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (44.4589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (46.7285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-18.9141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (44.8593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (181.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (123.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (146.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (383.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (392.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (280.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (65.9248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (175.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-10.5491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-3.18551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (51.8281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (159.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (459.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (388.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (404.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (383.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-58.3013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (200.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (440.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (390.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (540.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['O=C=O(1731)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['[CH2]C1OC1C([O])=O(14648)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['C=CC1OC1([O])[O](14609)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(138.019,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['[CH2]C1OC(=O)C1[O](14528)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.448e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C=CC(=O)C([O])=O(14649)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;HJ] for rate rule [CO-DeDe_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(64)', '[O]C(=O)C=O(3180)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0131003,'m^3/(mol*s)'), n=2.40999, Ea=(12.7705,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CdsJ-H] for rate rule [CO-DeH_O;CdsJ-H]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][C]=O(2059)', 'C=CC=O(5269)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.5e+06,'cm^3/(mol*s)'), n=2.16, Ea=(19.2464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdH_O;YJ] for rate rule [CO-CdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C=O(1731)', '[CH2]C=C[O](5266)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(246.209,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CJ] for rate rule [CO2;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 241.2 to 246.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['C=CC(O)=C([O])[O](14650)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['C=[C]C(O)C([O])=O(14651)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_Cs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['C=CC([O])=C([O])O(14652)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.99823e+07,'s^-1'), n=1.57622, Ea=(111.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]C([O])C(=O)O(14653)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC(O)C([O])=O(14654)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['[CH]=CC([O])C(=O)O(14655)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.936e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(64)', '[O]C=C([O])[O](3183)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.88428e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C]=O(2059)', '[CH2]C=C[O](5266)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', 'C=CC([O])=C([O])[O](14656)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', 'C=[C]C([O])C([O])=O(14657)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH]=CC([O])C([O])=O(14658)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['C=CC([O])[C]1OO1(14659)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.29014e+11,'s^-1'), n=0.514092, Ea=(347.198,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['[O]C(=O)C1[CH]CO1(14660)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.15968e+08,'s^-1'), n=1.10215, Ea=(132.51,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['[CH2][CH]C1OOC1=O(14553)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.006e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['[O]C1[CH]COC1=O(14569)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.14356e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['C=CC(=O)C(=O)O(14661)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CC([O])=C([O])[O](12041)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C=C([O])C([O])[O](14662)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]C([O])C([O])[O](14663)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C([O])[C]([O])O(14664)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=CC([O])C([O])[O](14665)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=CC([O])[C]([O])O(14666)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=CC([O])C([O])=O(12739)'],
    products = ['C=CC1OOC1=O(14551)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=COC([O])=C[O](14667)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O][C]1OC[CH]C1[O](14668)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(T)(63)', 'C=C[CH]C([O])=O(5461)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(T)(63)', 'C=CC([O])[C]=O(14526)'],
    products = ['C=CC([O])C([O])=O(12739)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

network(
    label = '3498',
    isomers = [
        'C=CC([O])C([O])=O(12739)',
    ],
    reactants = [
        ('O=C=O(1731)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3498',
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

