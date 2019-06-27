species(
    label = '[O]C=CCC[C]([O])[O](12716)',
    structure = SMILES('[O]C=CCC[C]([O])[O]'),
    E0 = (188.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,360,370,350,436.918,436.918,436.919,436.919,1991.29,1991.29],'cm^-1')),
        HinderedRotor(inertia=(0.18398,'amu*angstrom^2'), symmetry=1, barrier=(24.9227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.067755,'amu*angstrom^2'), symmetry=1, barrier=(9.17841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0677553,'amu*angstrom^2'), symmetry=1, barrier=(9.17835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14342,0.0686964,-7.24885e-05,4.54419e-08,-1.23597e-11,22825.2,31.4463], Tmin=(100,'K'), Tmax=(864.372,'K')), NASAPolynomial(coeffs=[7.89912,0.0374331,-1.82346e-05,3.59673e-09,-2.5677e-13,21657.3,-0.161317], Tmin=(864.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(C=COJ)"""),
)

species(
    label = 'C=C([O])[O](1172)',
    structure = SMILES('C=C([O])[O]'),
    E0 = (-40.8548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,219.703,219.72],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3790.78,'J/mol'), sigma=(6.02099,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.11 K, Pc=39.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85742,0.0201906,-6.55878e-06,-8.62627e-09,5.39429e-12,-4868.13,11.9276], Tmin=(100,'K'), Tmax=(974.512,'K')), NASAPolynomial(coeffs=[9.20068,0.00568109,-1.96826e-06,3.71357e-10,-2.78204e-14,-6651.8,-21.3196], Tmin=(974.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.8548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[O][CH]C1CCC1([O])[O](16918)',
    structure = SMILES('[O][CH]C1CCC1([O])[O]'),
    E0 = (306.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48259,0.0489927,-2.5598e-05,5.57596e-09,-4.55388e-13,36909.4,23.6408], Tmin=(100,'K'), Tmax=(2850.43,'K')), NASAPolynomial(coeffs=[26.8856,0.014748,-7.57714e-06,1.36119e-09,-8.57268e-14,22997.6,-119.651], Tmin=(2850.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclobutane) + radical(CCsJOH) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ) + radical(CCOJ)"""),
)

species(
    label = '[O][CH]C1CC[C]([O])O1(16919)',
    structure = SMILES('[O][CH]C1CC[C]([O])O1'),
    E0 = (241.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45343,0.0573952,-4.96223e-05,2.88354e-08,-7.52312e-12,29185.6,26.3831], Tmin=(100,'K'), Tmax=(894.201,'K')), NASAPolynomial(coeffs=[6.15573,0.0363609,-1.43384e-05,2.53024e-09,-1.68868e-13,28344.6,4.22293], Tmin=(894.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Tetrahydrofuran) + radical(CCsJOH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O]C=CCC=C([O])[O](16920)',
    structure = SMILES('[O]C=CCC=C([O])[O]'),
    E0 = (-61.5844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,268.981,268.996,269.002,269.004,269.013],'cm^-1')),
        HinderedRotor(inertia=(0.431275,'amu*angstrom^2'), symmetry=1, barrier=(22.1463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.43129,'amu*angstrom^2'), symmetry=1, barrier=(22.1471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.485949,0.0646582,-3.92674e-05,-6.38644e-09,9.9385e-12,-7268.97,29.1724], Tmin=(100,'K'), Tmax=(965.519,'K')), NASAPolynomial(coeffs=[20.1342,0.0126036,-3.98651e-06,7.31362e-10,-5.51373e-14,-12430.9,-72.0127], Tmin=(965.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.5844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C]([O])CCC=C=O(16921)',
    structure = SMILES('[O][C]([O])CCC=C=O'),
    E0 = (171.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2120,512.5,787.5,3010,987.5,1337.5,450,1655,360,370,350,180,180,1597.45,1601.81,1606.15],'cm^-1')),
        HinderedRotor(inertia=(0.217325,'amu*angstrom^2'), symmetry=1, barrier=(4.99672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215624,'amu*angstrom^2'), symmetry=1, barrier=(4.95761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216864,'amu*angstrom^2'), symmetry=1, barrier=(4.98613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0614,0.0771198,-0.000132443,1.33671e-07,-5.17297e-11,20676.4,31.1801], Tmin=(100,'K'), Tmax=(828.473,'K')), NASAPolynomial(coeffs=[1.82865,0.0451055,-2.3223e-05,4.53619e-09,-3.15726e-13,21520.8,33.4864], Tmin=(828.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-(Cdd-O2d)CsH) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P)"""),
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
    label = '[O]C=CC[CH]C([O])[O](16922)',
    structure = SMILES('[O]C=CC[CH]C([O])[O]'),
    E0 = (183.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33451,0.0625722,-5.37218e-05,2.51648e-08,-5.05496e-12,22177.2,32.5533], Tmin=(100,'K'), Tmax=(1139.73,'K')), NASAPolynomial(coeffs=[9.48803,0.0339561,-1.60596e-05,3.1345e-09,-2.22539e-13,20318.7,-7.84895], Tmin=(1139.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCOJ) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = '[O][C]([O])CCC=[C]O(16923)',
    structure = SMILES('[O][C]([O])CCC=[C]O'),
    E0 = (287.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3615,1277.5,1000,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596776,0.081827,-0.000119857,1.02852e-07,-3.57836e-11,34663.8,35.105], Tmin=(100,'K'), Tmax=(798.116,'K')), NASAPolynomial(coeffs=[7.89958,0.0363038,-1.75287e-05,3.36978e-09,-2.34048e-13,33782.2,3.30037], Tmin=(798.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C=C[CH]CC([O])[O](16924)',
    structure = SMILES('[O]C=C[CH]CC([O])[O]'),
    E0 = (124.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40757,0.0617533,-5.05964e-05,2.25352e-08,-4.35878e-12,15103.1,28.6871], Tmin=(100,'K'), Tmax=(1164.15,'K')), NASAPolynomial(coeffs=[8.97741,0.0357439,-1.7084e-05,3.34418e-09,-2.37594e-13,13340.6,-8.98366], Tmin=(1164.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(CCOJ) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = '[O]C=CC[CH][C]([O])O(16925)',
    structure = SMILES('[O]C=CC[CH][C]([O])O'),
    E0 = (163.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646086,0.0704065,-6.74314e-05,3.30019e-08,-6.45491e-12,19747.5,35.5394], Tmin=(100,'K'), Tmax=(1230.42,'K')), NASAPolynomial(coeffs=[15.14,0.023288,-9.9894e-06,1.87869e-09,-1.31204e-13,16180.8,-37.3906], Tmin=(1230.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(Cs_P) + radical(CCJCO) + radical(C=COJ)"""),
)

species(
    label = '[O][C]([O])CC[C]=CO(16926)',
    structure = SMILES('[O][C]([O])CC[C]=CO'),
    E0 = (285.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3615,1277.5,1000,360,370,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.474661,0.0830687,-0.000111241,8.33361e-08,-2.54528e-11,34441,32.8995], Tmin=(100,'K'), Tmax=(796.183,'K')), NASAPolynomial(coeffs=[10.6196,0.0321,-1.52156e-05,2.92982e-09,-2.05e-13,32825.6,-13.7316], Tmin=(796.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C=[C]CCC([O])[O](16927)',
    structure = SMILES('[O]C=[C]CCC([O])[O]'),
    E0 = (221.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,385.239,385.239,385.244,385.245,2128.37,2128.37],'cm^-1')),
        HinderedRotor(inertia=(0.259905,'amu*angstrom^2'), symmetry=1, barrier=(27.3718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0971415,'amu*angstrom^2'), symmetry=1, barrier=(10.2308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0971433,'amu*angstrom^2'), symmetry=1, barrier=(10.2307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26351,0.0690599,-6.61916e-05,1.6321e-08,1.51767e-11,26737.6,30.8674], Tmin=(100,'K'), Tmax=(537.678,'K')), NASAPolynomial(coeffs=[7.0885,0.0389595,-1.91385e-05,3.75715e-09,-2.66278e-13,25919.9,4.60059], Tmin=(537.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C[CH]C[C]([O])O(16928)',
    structure = SMILES('[O]C=C[CH]C[C]([O])O'),
    E0 = (104.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.709632,0.0696875,-6.46013e-05,3.06896e-08,-5.87017e-12,12673.9,31.7082], Tmin=(100,'K'), Tmax=(1249.29,'K')), NASAPolynomial(coeffs=[14.6525,0.0250448,-1.09993e-05,2.08553e-09,-1.4606e-13,9190.17,-38.6612], Tmin=(1249.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(CCOJ) + radical(Cs_P) + radical(C=COJ)"""),
)

species(
    label = '[O][C]([O])C[CH]C=CO(16929)',
    structure = SMILES('[O][C]([O])C[CH]C=CO'),
    E0 = (188.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82601,0.0718738,-7.28412e-05,3.94997e-08,-8.79072e-12,22797.9,30.0546], Tmin=(100,'K'), Tmax=(1071.41,'K')), NASAPolynomial(coeffs=[12.1661,0.0295373,-1.35695e-05,2.61918e-09,-1.85184e-13,20367.9,-25.4368], Tmin=(1071.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_P) + radical(CCOJ) + radical(Allyl_S) + radical(CCOJ)"""),
)

species(
    label = '[O]C([O])CC[CH][C]=O(16930)',
    structure = SMILES('[O]C([O])CC[CH][C]=O'),
    E0 = (167.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.899193,0.0797546,-0.000130048,1.26585e-07,-4.77046e-11,20301.2,32.4778], Tmin=(100,'K'), Tmax=(837.736,'K')), NASAPolynomial(coeffs=[2.75397,0.0452393,-2.23031e-05,4.28015e-09,-2.95036e-13,20890.8,29.2319], Tmin=(837.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCOJ) + radical(CCJCHO) + radical(CCOJ)"""),
)

species(
    label = '[O]C=[C]CC[C]([O])O(16931)',
    structure = SMILES('[O]C=[C]CC[C]([O])O'),
    E0 = (201.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3615,1277.5,1000,360,370,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.670224,0.0770926,-8.92774e-05,5.60923e-08,-1.43776e-11,24303.7,33.4425], Tmin=(100,'K'), Tmax=(941.028,'K')), NASAPolynomial(coeffs=[11.7948,0.0298048,-1.38991e-05,2.68997e-09,-1.90097e-13,22210.1,-19.5508], Tmin=(941.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O][C]([O])[CH]CC=CO(16932)',
    structure = SMILES('[O][C]([O])[CH]CC=CO'),
    E0 = (247.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,360,370,350,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.762344,0.0726064,-7.576e-05,4.19608e-08,-9.445e-12,29871.5,33.8852], Tmin=(100,'K'), Tmax=(1066.46,'K')), NASAPolynomial(coeffs=[12.7621,0.0275987,-1.24557e-05,2.38795e-09,-1.68313e-13,27312,-24.7786], Tmin=(1066.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = '[O][C](O)CC[CH][C]=O(16933)',
    structure = SMILES('[O][C](O)CC[CH][C]=O'),
    E0 = (147.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.595891,0.0829478,-0.000127027,1.11902e-07,-3.90669e-11,17854.9,34.0889], Tmin=(100,'K'), Tmax=(833.688,'K')), NASAPolynomial(coeffs=[7.24292,0.0365222,-1.73468e-05,3.28546e-09,-2.25235e-13,17251.6,6.25899], Tmin=(833.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Cs_P) + radical(CCCJ=O) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]([O])[O](2304)',
    structure = SMILES('[CH2][C]([O])[O]'),
    E0 = (411.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,1740.05,1740.26],'cm^-1')),
        HinderedRotor(inertia=(0.00349141,'amu*angstrom^2'), symmetry=1, barrier=(7.50228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01691,0.0315826,-7.13533e-05,8.28598e-08,-3.37644e-11,49511,16.8536], Tmin=(100,'K'), Tmax=(857.006,'K')), NASAPolynomial(coeffs=[-0.765291,0.0236786,-1.27872e-05,2.50409e-09,-1.7282e-13,51097.8,39.9925], Tmin=(857.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCO) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[CH2]C[C]([O])[O](2320)',
    structure = SMILES('[CH2]C[C]([O])[O]'),
    E0 = (381.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,360,370,350,1619.61,1622.32,1624.36],'cm^-1')),
        HinderedRotor(inertia=(0.297191,'amu*angstrom^2'), symmetry=1, barrier=(6.833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00366612,'amu*angstrom^2'), symmetry=1, barrier=(6.85169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5224,0.0423956,-7.78802e-05,8.71711e-08,-3.58832e-11,45905.9,22.0246], Tmin=(100,'K'), Tmax=(836.967,'K')), NASAPolynomial(coeffs=[-0.833248,0.0335228,-1.73354e-05,3.38608e-09,-2.35344e-13,47340.1,42.8286], Tmin=(836.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ) + radical(RCCJ)"""),
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
    label = '[O]C=C[CH]C[C]([O])[O](16934)',
    structure = SMILES('[O]C=C[CH]C[C]([O])[O]'),
    E0 = (330.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,519.35,519.445,519.481,519.492,519.541,519.844],'cm^-1')),
        HinderedRotor(inertia=(0.00904029,'amu*angstrom^2'), symmetry=1, barrier=(1.73176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00904144,'amu*angstrom^2'), symmetry=1, barrier=(1.73221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00903449,'amu*angstrom^2'), symmetry=1, barrier=(1.73197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39115,0.0629275,-6.18079e-05,3.4913e-08,-8.57506e-12,39788.6,29.448], Tmin=(100,'K'), Tmax=(945.928,'K')), NASAPolynomial(coeffs=[8.06192,0.0347187,-1.70751e-05,3.38597e-09,-2.42599e-13,38526.6,-2.36368], Tmin=(945.928,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(CCOJ) + radical(Allyl_S) + radical(Cs_P)"""),
)

species(
    label = '[O]C=CC[CH][C]([O])[O](16935)',
    structure = SMILES('[O]C=CC[CH][C]([O])[O]'),
    E0 = (388.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,354.944,354.946,354.985,724.922,1904.66,1904.79],'cm^-1')),
        HinderedRotor(inertia=(0.00133818,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00133761,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25588,'amu*angstrom^2'), symmetry=1, barrier=(22.8727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33902,0.0635399,-6.43734e-05,3.69923e-08,-9.09161e-12,46861.7,33.2361], Tmin=(100,'K'), Tmax=(954.381,'K')), NASAPolynomial(coeffs=[8.69437,0.0327112,-1.59186e-05,3.1441e-09,-2.2481e-13,45457.7,-1.90559], Tmin=(954.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCJCO) + radical(CCOJ) + radical(C=COJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C=[C]CC[C]([O])[O](16936)',
    structure = SMILES('[O]C=[C]CC[C]([O])[O]'),
    E0 = (426.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,360,370,350,395.564,395.564,395.565,395.565,2063.09,2063.09],'cm^-1')),
        HinderedRotor(inertia=(0.0615264,'amu*angstrom^2'), symmetry=1, barrier=(6.8316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0615263,'amu*angstrom^2'), symmetry=1, barrier=(6.8316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239962,'amu*angstrom^2'), symmetry=1, barrier=(26.6443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56793,0.0647495,-4.48585e-05,-5.05981e-08,7.73238e-11,51409,30.5584], Tmin=(100,'K'), Tmax=(471.908,'K')), NASAPolynomial(coeffs=[7.42719,0.0356528,-1.77485e-05,3.46088e-09,-2.42485e-13,50627,4.26445], Tmin=(471.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCOJ) + radical(Cds_S) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O][C]([O])CC[CH][C]=O(16937)',
    structure = SMILES('[O][C]([O])CC[CH][C]=O'),
    E0 = (373.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,1855,455,950,3025,407.5,1350,352.5,253.027,866.068,1154.76,1510.09,1881.89],'cm^-1')),
        HinderedRotor(inertia=(0.102891,'amu*angstrom^2'), symmetry=1, barrier=(3.27266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102891,'amu*angstrom^2'), symmetry=1, barrier=(3.27266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102891,'amu*angstrom^2'), symmetry=1, barrier=(3.27266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102891,'amu*angstrom^2'), symmetry=1, barrier=(3.27266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924411,0.0806277,-0.00014096,1.39349e-07,-5.22965e-11,44984.3,33.0752], Tmin=(100,'K'), Tmax=(848.838,'K')), NASAPolynomial(coeffs=[2.62266,0.0428136,-2.14573e-05,4.11894e-09,-2.82866e-13,45770,31.4867], Tmin=(848.838,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCCJ=O) + radical(Cs_P) + radical(CCOJ) + radical(CCJCHO)"""),
)

species(
    label = '[O][C]([O])CCC1[CH]O1(16938)',
    structure = SMILES('[O][C]([O])CCC1[CH]O1'),
    E0 = (329.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635587,0.0786552,-0.000105698,8.32563e-08,-2.60605e-11,39725.9,31.1056], Tmin=(100,'K'), Tmax=(906.702,'K')), NASAPolynomial(coeffs=[8.45259,0.0334753,-1.32628e-05,2.28318e-09,-1.47475e-13,38748,-3.4171], Tmin=(906.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Ethylene_oxide) + radical(Cs_P) + radical(CCsJO) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C1[CH]CCC1([O])[O](16939)',
    structure = SMILES('[O]C1[CH]CCC1([O])[O]'),
    E0 = (245.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72202,0.0378874,1.68078e-05,-4.28714e-08,1.7134e-11,29614.8,28.4826], Tmin=(100,'K'), Tmax=(1036.69,'K')), NASAPolynomial(coeffs=[11.2973,0.0294562,-1.22507e-05,2.34703e-09,-1.68666e-13,26097.2,-25.4477], Tmin=(1036.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CC(C)(O)OJ) + radical(CCJCO) + radical(CC(C)(O)OJ)"""),
)

species(
    label = '[O][C]1CC[CH]C([O])O1(16859)',
    structure = SMILES('[O][C]1CC[CH]C([O])O1'),
    E0 = (211.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8291,0.0408221,-1.32901e-06,-1.85135e-08,8.01552e-12,25573.8,27.8307], Tmin=(100,'K'), Tmax=(1043.97,'K')), NASAPolynomial(coeffs=[7.23784,0.0345791,-1.31651e-05,2.33151e-09,-1.57979e-13,23655.4,-2.27543], Tmin=(1043.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(Oxane) + radical(CCOJ) + radical(Cs_P) + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = '[O]C([O])=CCC=CO(16940)',
    structure = SMILES('[O]C([O])=CCC=CO'),
    E0 = (-203.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0796723,0.0715672,-4.24582e-05,-1.28792e-08,1.48026e-11,-24266.3,29.2206], Tmin=(100,'K'), Tmax=(934.839,'K')), NASAPolynomial(coeffs=[23.3614,0.00896009,-1.38783e-06,1.82394e-10,-1.60168e-14,-30236.4,-90.1813], Tmin=(934.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C([O])CCC=C=O(16941)',
    structure = SMILES('[O]C([O])CCC=C=O'),
    E0 = (-34.1119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04019,0.0762012,-0.000121391,1.20779e-07,-4.71243e-11,-4006.95,30.5682], Tmin=(100,'K'), Tmax=(815.122,'K')), NASAPolynomial(coeffs=[1.92909,0.0475857,-2.41011e-05,4.70518e-09,-3.28553e-13,-3346.13,31.4039], Tmin=(815.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.1119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-(Cdd-O2d)CsH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O][CH]C[CH]C[C]([O])[O](16942)',
    structure = SMILES('[O][CH]C[CH]C[C]([O])[O]'),
    E0 = (569.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3050,390,425,1340,1360,335,370,191.444,273.94,842.63,1036.63,2170.29,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0457736,'amu*angstrom^2'), symmetry=1, barrier=(1.0669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0457736,'amu*angstrom^2'), symmetry=1, barrier=(1.0669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0457736,'amu*angstrom^2'), symmetry=1, barrier=(1.0669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0457736,'amu*angstrom^2'), symmetry=1, barrier=(1.0669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738973,0.098358,-0.000211395,2.32876e-07,-9.142e-11,68544.5,36.8047], Tmin=(100,'K'), Tmax=(867.302,'K')), NASAPolynomial(coeffs=[-5.86844,0.0618493,-3.24074e-05,6.24751e-09,-4.26388e-13,72209.9,82.2644], Tmin=(867.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(RCCJCC) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O][CH]CC[CH][C]([O])[O](16943)',
    structure = SMILES('[O][CH]CC[CH][C]([O])[O]'),
    E0 = (574.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3050,390,425,1340,1360,335,370,197.93,1210.3,1219.81,1269.84,1451.25,1647.64,1708.48],'cm^-1')),
        HinderedRotor(inertia=(0.11364,'amu*angstrom^2'), symmetry=1, barrier=(3.04142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11364,'amu*angstrom^2'), symmetry=1, barrier=(3.04142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11364,'amu*angstrom^2'), symmetry=1, barrier=(3.04142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11364,'amu*angstrom^2'), symmetry=1, barrier=(3.04142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.656573,0.0929492,-0.000179856,1.8869e-07,-7.29539e-11,69209.1,37.0329], Tmin=(100,'K'), Tmax=(853.731,'K')), NASAPolynomial(coeffs=[-0.960302,0.0534663,-2.78037e-05,5.38891e-09,-3.70812e-13,71200.1,54.6215], Tmin=(853.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH) + radical(CCJCO) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]C[CH]C[CH][C]([O])[O](16944)',
    structure = SMILES('[O]C[CH]C[CH][C]([O])[O]'),
    E0 = (594.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3050,390,425,1340,1360,335,370,189.082,405.804,808.226,1837.37,2582.07,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0114318,'amu*angstrom^2'), symmetry=1, barrier=(0.267451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114318,'amu*angstrom^2'), symmetry=1, barrier=(0.267451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114318,'amu*angstrom^2'), symmetry=1, barrier=(0.267451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114318,'amu*angstrom^2'), symmetry=1, barrier=(0.267451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08765,0.0794027,-0.000141579,1.492e-07,-5.90666e-11,71555.2,38.4321], Tmin=(100,'K'), Tmax=(837.204,'K')), NASAPolynomial(coeffs=[-0.949901,0.0522563,-2.6862e-05,5.23136e-09,-3.62943e-13,73188.9,55.6193], Tmin=(837.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(594.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCJCO) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C([O])([O])CC=C[O](12712)',
    structure = SMILES('[CH2]C([O])([O])CC=C[O]'),
    E0 = (188.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,180,180,701.297,1516.59,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149164,'amu*angstrom^2'), symmetry=1, barrier=(3.42958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149164,'amu*angstrom^2'), symmetry=1, barrier=(3.42958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149164,'amu*angstrom^2'), symmetry=1, barrier=(3.42958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.829728,0.0754463,-8.4463e-05,5.40986e-08,-1.4621e-11,22731.7,29.2922], Tmin=(100,'K'), Tmax=(880.868,'K')), NASAPolynomial(coeffs=[9.48471,0.0361442,-1.75369e-05,3.44698e-09,-2.45463e-13,21206.9,-11.3652], Tmin=(880.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)(O)OJ) + radical(CJC(O)2C) + radical(CC(C)(O)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C1([O])CCC=CO1(16945)',
    structure = SMILES('[O]C1([O])CCC=CO1'),
    E0 = (-83.5851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71477,0.040608,6.26938e-06,-3.33691e-08,1.50006e-11,-9962.1,22.4398], Tmin=(100,'K'), Tmax=(976.175,'K')), NASAPolynomial(coeffs=[10.1877,0.0286393,-1.02976e-05,1.81928e-09,-1.25485e-13,-12700.3,-23.7848], Tmin=(976.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.5851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C=CCCC([O])=O(12726)',
    structure = SMILES('[O]C=CCCC([O])=O'),
    E0 = (-230.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22315,0.0601708,-4.35576e-05,1.58483e-08,-2.38339e-12,-27612.1,27.0637], Tmin=(100,'K'), Tmax=(1506.14,'K')), NASAPolynomial(coeffs=[12.5793,0.030011,-1.35208e-05,2.55295e-09,-1.76536e-13,-31032.9,-32.3742], Tmin=(1506.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ)"""),
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
    label = '[O][C][O](2319)',
    structure = SMILES('[O][C][O]'),
    E0 = (504.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([786.798,787.166,787.343],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.36312,0.000897068,5.83704e-06,-5.0699e-09,1.05455e-12,60693,5.58445], Tmin=(100,'K'), Tmax=(1870.77,'K')), NASAPolynomial(coeffs=[6.74008,0.00434757,-3.77128e-06,7.92205e-10,-5.4643e-14,58310.5,-11.3625], Tmin=(1870.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(CH2_triplet) + radical(OCOJ) + radical(OCOJ)"""),
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
    label = '[CH]=CCC[C]([O])[O](16946)',
    structure = SMILES('[CH]=CCC[C]([O])[O]'),
    E0 = (503.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,360,370,350,180,180,1804.66,1804.77],'cm^-1')),
        HinderedRotor(inertia=(0.12521,'amu*angstrom^2'), symmetry=1, barrier=(2.87883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125192,'amu*angstrom^2'), symmetry=1, barrier=(2.87841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125174,'amu*angstrom^2'), symmetry=1, barrier=(2.878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34029,0.0680222,-0.000106209,1.04772e-07,-4.08043e-11,60630.6,29.8471], Tmin=(100,'K'), Tmax=(812.004,'K')), NASAPolynomial(coeffs=[2.35861,0.0425699,-2.14408e-05,4.18218e-09,-2.92123e-13,61138.9,29.2949], Tmin=(812.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(503.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O][C]([O])CC=CC=O(16947)',
    structure = SMILES('[O][C]([O])CC=CC=O'),
    E0 = (156.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,360,370,350,180,180,1438.59,4000],'cm^-1')),
        HinderedRotor(inertia=(0.185209,'amu*angstrom^2'), symmetry=1, barrier=(4.25831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.876045,'amu*angstrom^2'), symmetry=1, barrier=(20.142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.013641,'amu*angstrom^2'), symmetry=1, barrier=(20.1497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.40718,0.0408889,-2.04174e-05,3.34452e-09,-1.45372e-13,18765.9,18.9699], Tmin=(100,'K'), Tmax=(2888.71,'K')), NASAPolynomial(coeffs=[63.3536,-0.025101,5.73095e-06,-8.16634e-10,5.26155e-14,-21812.6,-339.233], Tmin=(2888.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O][C]([O])C[CH]CC=O(16885)',
    structure = SMILES('[O][C]([O])C[CH]CC=O'),
    E0 = (245.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,220.295,636.108,807.439,1925.79,3696.57],'cm^-1')),
        HinderedRotor(inertia=(0.0144873,'amu*angstrom^2'), symmetry=1, barrier=(0.447364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144873,'amu*angstrom^2'), symmetry=1, barrier=(0.447364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144873,'amu*angstrom^2'), symmetry=1, barrier=(0.447364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144873,'amu*angstrom^2'), symmetry=1, barrier=(0.447364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03732,0.07874,-0.000134783,1.38705e-07,-5.45e-11,29634.3,34.5843], Tmin=(100,'K'), Tmax=(829.474,'K')), NASAPolynomial(coeffs=[0.442209,0.0499874,-2.56022e-05,4.99361e-09,-3.47379e-13,30820.9,43.9017], Tmin=(829.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P) + radical(CCJCC=O)"""),
)

species(
    label = '[O][C]([O])CCC[C]=O(12229)',
    structure = SMILES('[O][C]([O])CCC[C]=O'),
    E0 = (205.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1855,455,950,360,370,350,180,180,375.395,1838.84,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0429464,'amu*angstrom^2'), symmetry=1, barrier=(0.987423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0429464,'amu*angstrom^2'), symmetry=1, barrier=(0.987423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0429464,'amu*angstrom^2'), symmetry=1, barrier=(0.987423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0429464,'amu*angstrom^2'), symmetry=1, barrier=(0.987423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772279,0.0868625,-0.000158062,1.61535e-07,-6.19134e-11,24837.9,33.8592], Tmin=(100,'K'), Tmax=(847.833,'K')), NASAPolynomial(coeffs=[0.801269,0.0492317,-2.515e-05,4.86284e-09,-3.35088e-13,26180.5,41.6713], Tmin=(847.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O][C]([O])[CH]CCC=O(16948)',
    structure = SMILES('[O][C]([O])[CH]CCC=O'),
    E0 = (245.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,360,370,350,3025,407.5,1350,352.5,220.295,636.108,807.439,1925.79,3696.57],'cm^-1')),
        HinderedRotor(inertia=(0.0144873,'amu*angstrom^2'), symmetry=1, barrier=(0.447364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144873,'amu*angstrom^2'), symmetry=1, barrier=(0.447364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144873,'amu*angstrom^2'), symmetry=1, barrier=(0.447364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0144873,'amu*angstrom^2'), symmetry=1, barrier=(0.447364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03732,0.07874,-0.000134783,1.38705e-07,-5.45e-11,29634.3,34.5843], Tmin=(100,'K'), Tmax=(829.474,'K')), NASAPolynomial(coeffs=[0.442209,0.0499874,-2.56022e-05,4.99361e-09,-3.47379e-13,30820.9,43.9017], Tmin=(829.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(Cs_P) + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = '[O]C1([O])CC[CH][CH]O1(16949)',
    structure = SMILES('[O]C1([O])CC[CH][CH]O1'),
    E0 = (206.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68425,0.048258,-2.04459e-05,-4.35441e-09,5.52704e-12,24981.9,25.8095], Tmin=(100,'K'), Tmax=(860.64,'K')), NASAPolynomial(coeffs=[7.00432,0.0333144,-1.14506e-05,1.88472e-09,-1.21711e-13,23703.9,-1.16314], Tmin=(860.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Oxane) + radical(CCOJ) + radical(CCsJOCs) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = '[O][C]1CC[CH][CH]OO1(16795)',
    structure = SMILES('[O][C]1CC[CH][CH]OO1'),
    E0 = (451.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18216,0.0216061,0.000131891,-2.15996e-07,9.29006e-11,54414.6,31.5897], Tmin=(100,'K'), Tmax=(914.768,'K')), NASAPolynomial(coeffs=[32.3678,-0.00965786,1.0816e-05,-2.1592e-09,1.35263e-13,44311.6,-140.121], Tmin=(914.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(Cs_P) + radical(CCOJ) + radical(CCsJOOC)"""),
)

species(
    label = '[O]C([O])CC=CC=O(16950)',
    structure = SMILES('[O]C([O])CC=CC=O'),
    E0 = (-48.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95153,0.0443128,-2.1736e-05,3.54668e-09,-1.47942e-13,-5895.3,19.9819], Tmin=(100,'K'), Tmax=(2773.57,'K')), NASAPolynomial(coeffs=[54.7306,-0.0140493,1.78505e-06,-1.73834e-10,1.31706e-14,-39782.9,-287.11], Tmin=(2773.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C([O])=CCCC=O(16951)',
    structure = SMILES('[O]C([O])=CCCC=O'),
    E0 = (-206.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0337,0.0693914,-7.03494e-05,3.98451e-08,-9.49822e-12,-24685.4,27.4823], Tmin=(100,'K'), Tmax=(991.151,'K')), NASAPolynomial(coeffs=[9.99142,0.0332405,-1.56388e-05,3.04555e-09,-2.16166e-13,-26461.1,-15.6538], Tmin=(991.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C=O)C[C]([O])[O](12718)',
    structure = SMILES('[CH2]C(C=O)C[C]([O])[O]'),
    E0 = (250.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,360,370,350,3000,3100,440,815,1455,1000,180,1640.12,1640.49,1640.57],'cm^-1')),
        HinderedRotor(inertia=(0.224676,'amu*angstrom^2'), symmetry=1, barrier=(5.16574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224742,'amu*angstrom^2'), symmetry=1, barrier=(5.16726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2247,'amu*angstrom^2'), symmetry=1, barrier=(5.1663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224705,'amu*angstrom^2'), symmetry=1, barrier=(5.1664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4578.58,'J/mol'), sigma=(7.43252,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=715.16 K, Pc=25.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.566852,0.0918701,-0.000168678,1.70939e-07,-6.50303e-11,30231.8,33.991], Tmin=(100,'K'), Tmax=(847.468,'K')), NASAPolynomial(coeffs=[1.51395,0.0491426,-2.53367e-05,4.9111e-09,-3.38681e-13,31445.1,37.684], Tmin=(847.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(Cs_P) + radical(CJC(C)C=O) + radical(CCOJ)"""),
)

species(
    label = '[O]C1([O])CCC1C=O(12720)',
    structure = SMILES('[O]C1([O])CCC1C=O'),
    E0 = (-19.1241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5601,0.0490294,-2.29311e-05,2.75923e-09,3.80376e-13,-2208.48,26.136], Tmin=(100,'K'), Tmax=(1597.16,'K')), NASAPolynomial(coeffs=[13.7296,0.0283658,-1.27415e-05,2.35334e-09,-1.58304e-13,-7347.53,-42.1914], Tmin=(1597.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.1241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
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
    label = '[CH]CC[C]([O])[O](2321)',
    structure = SMILES('[CH]CC[C]([O])[O]'),
    E0 = (600.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,208.297,987.959,1065.73,1204.48,1376.55,1609.16,1857.45],'cm^-1')),
        HinderedRotor(inertia=(0.124104,'amu*angstrom^2'), symmetry=1, barrier=(3.30006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124104,'amu*angstrom^2'), symmetry=1, barrier=(3.30006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124104,'amu*angstrom^2'), symmetry=1, barrier=(3.30006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73906,0.0600761,-0.000102933,1.05719e-07,-4.14362e-11,72295.9,26.2983], Tmin=(100,'K'), Tmax=(830.957,'K')), NASAPolynomial(coeffs=[1.33876,0.0378971,-1.9382e-05,3.77563e-09,-2.62397e-13,73194.7,33.1632], Tmin=(830.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCOJ) + radical(CCJ2_triplet) + radical(Cs_P)"""),
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
    E0 = (188.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (311.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (292.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (188.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (385.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (188.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (345.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (450.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (313.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (344.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (435.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (322.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (268.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (302.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (300.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (352.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (342.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (261.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (501.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (603.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (547.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (600.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (639.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (585.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (402.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (247.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (243.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (213.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (213.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (591.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (637.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (619.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (345.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (196.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (188.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (660.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (910.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (368.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (334.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (370.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (364.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (370.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (259.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (451.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (267.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (277.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (410.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (196.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (668.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['C=C([O])[O](1172)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O][CH]C1CCC1([O])[O](16918)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O][CH]C1CC[C]([O])O1(16919)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.84666e+10,'s^-1'), n=0.385799, Ea=(103.916,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R6;doublebond_intra;radadd_intra] for rate rule [R6;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[O]C=CCC=C([O])[O](16920)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(38.7458,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 33.7 to 38.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[O][C]([O])CCC=C=O(16921)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])[O](1172)', '[CH2]C=C[O](5266)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.000178913,'m^3/(mol*s)'), n=2.74787, Ea=(139.528,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-CdHH]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 135.1 to 139.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C=CC[CH]C([O])[O](16922)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C]([O])CCC=[C]O(16923)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C=C[CH]CC([O])[O](16924)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C=CC[CH][C]([O])O(16925)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.61628e+09,'s^-1'), n=1.19923, Ea=(155.469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]([O])CC[C]=CO(16926)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C=[C]CCC([O])[O](16927)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C=C[CH]C[C]([O])O(16928)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.21784e+06,'s^-1'), n=1.98643, Ea=(79.1566,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4HJ_1;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O][C]([O])C[CH]C=CO(16929)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 289 used for R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C([O])CC[CH][C]=O(16930)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.234e+06,'s^-1'), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SSSD;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C=[C]CC[C]([O])O(16931)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C]([O])[CH]CC=CO(16932)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O][C](O)CC[CH][C]=O(16933)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.48874e+06,'s^-1'), n=0.972854, Ea=(72.9565,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;O_rad_out;XH_out] for rate rule [R6HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][C]([O])[O](2304)', '[CH2]C=C[O](5266)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C[C]([O])[O](2320)', '[CH]=C[O](602)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.00218e+08,'m^3/(mol*s)'), n=-0.446058, Ea=(0.74957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H2/Cs] + [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/Cs]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[O]C=C[CH]C[C]([O])[O](16934)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.66284e+07,'m^3/(mol*s)'), n=0.108445, Ea=(5.74998,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[O]C=CC[CH][C]([O])[O](16935)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[O]C=[C]CC[C]([O])[O](16936)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[O][C]([O])CC[CH][C]=O(16937)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O][C]([O])CCC1[CH]O1(16938)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C1[CH]CCC1([O])[O](16939)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.88298e+06,'s^-1'), n=1.13333, Ea=(58.2971,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_CsCs_HH_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O][C]1CC[CH]C([O])O1(16859)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.37086e+10,'s^-1'), n=0.209443, Ea=(54.3406,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C([O])=CCC=CO(16940)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C([O])CCC=C=O(16941)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O][CH]C[CH]C[C]([O])[O](16942)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O][CH]CC[CH][C]([O])[O](16943)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C[CH]C[CH][C]([O])[O](16944)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])([O])CC=C[O](12712)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C1([O])CCC=CO1(16945)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C=CCCC([O])=O(12726)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]CC=C[O](743)', '[O][C][O](2319)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['O(T)(63)', '[CH]=CCC[C]([O])[O](16946)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(8)', '[O][C]([O])CC=CC=O(16947)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.76955,'m^3/(mol*s)'), n=1.94497, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-OneDeH;HJ] for rate rule [Cds-CsH_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][C]([O])[O](2304)', 'C=CC=O(5269)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.0102751,'m^3/(mol*s)'), n=2.40501, Ea=(4.48561,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDeH;CJ] for rate rule [Cds-HH_Cds-COH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O][C]([O])C[CH]CC=O(16885)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O][C]([O])CCC[C]=O(12229)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O][C]([O])[CH]CCC=O(16948)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C1([O])CC[CH][CH]O1(16949)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(8.58001e+07,'s^-1'), n=0.730566, Ea=(70.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O][C]1CC[CH][CH]OO1(16795)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.1298e+13,'s^-1'), n=0.287, Ea=(262.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 256.5 to 262.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C([O])CC=CC=O(16950)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C([O])=CCCC=O(16951)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=O)C[C]([O])[O](12718)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[O]C=CCC[C]([O])[O](12716)'],
    products = ['[O]C1([O])CCC1C=O(12720)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=O(373)', '[CH]CC[C]([O])[O](2321)'],
    products = ['[O]C=CCC[C]([O])[O](12716)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '3471',
    isomers = [
        '[O]C=CCC[C]([O])[O](12716)',
    ],
    reactants = [
        ('C=C([O])[O](1172)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3471',
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

