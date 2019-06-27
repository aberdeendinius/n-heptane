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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00748,0.0534506,-1.61976e-05,-2.04947e-08,1.22809e-11,52644,30.7954], Tmin=(100,'K'), Tmax=(999.958,'K')), NASAPolynomial(coeffs=[16.1987,0.020262,-7.78247e-06,1.47611e-09,-1.07528e-13,48227.1,-49.3871], Tmin=(999.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
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
    label = 'C=[C]C=C(O)[C]=C(28823)',
    structure = SMILES('C=[C]C=C(O)[C]=C'),
    E0 = (329.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.41514,'amu*angstrom^2'), symmetry=1, barrier=(32.5369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41334,'amu*angstrom^2'), symmetry=1, barrier=(32.4955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41332,'amu*angstrom^2'), symmetry=1, barrier=(32.4951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.290419,0.0830004,-9.76143e-05,5.59586e-08,-1.20009e-11,39829,24.4572], Tmin=(100,'K'), Tmax=(1300.02,'K')), NASAPolynomial(coeffs=[20.0388,0.0096718,-5.69708e-07,-1.84336e-10,2.20773e-14,35454.1,-75.4507], Tmin=(1300.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC(O)[CH][C]=C(24342)',
    structure = SMILES('C#CC(O)[CH][C]=C'),
    E0 = (411.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2175,525,750,770,3400,2100,217.462,217.569],'cm^-1')),
        HinderedRotor(inertia=(0.935215,'amu*angstrom^2'), symmetry=1, barrier=(31.4584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.935691,'amu*angstrom^2'), symmetry=1, barrier=(31.4771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937814,'amu*angstrom^2'), symmetry=1, barrier=(31.4383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926413,'amu*angstrom^2'), symmetry=1, barrier=(31.4407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520871,0.0701316,-6.91644e-05,3.27962e-08,-5.49108e-12,49635.1,26.0843], Tmin=(100,'K'), Tmax=(1008.58,'K')), NASAPolynomial(coeffs=[16.2508,0.017783,-6.23546e-06,1.06644e-09,-7.16132e-14,45951.7,-52.4686], Tmin=(1008.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=CC(O)[C]=C(27788)',
    structure = SMILES('[CH]=C=CC(O)[C]=C'),
    E0 = (443.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.694272,'amu*angstrom^2'), symmetry=1, barrier=(15.9627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694783,'amu*angstrom^2'), symmetry=1, barrier=(15.9744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694709,'amu*angstrom^2'), symmetry=1, barrier=(15.9727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06512,0.0640079,-6.62022e-05,3.69969e-08,-8.35491e-12,53484.7,28.3977], Tmin=(100,'K'), Tmax=(1069.83,'K')), NASAPolynomial(coeffs=[11.8516,0.0236781,-9.65573e-06,1.75967e-09,-1.20544e-13,51176.7,-24.3687], Tmin=(1069.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
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
    label = 'C=[C]C=CO(19110)',
    structure = SMILES('C=[C]C=CO'),
    E0 = (84.1241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.37885,'amu*angstrom^2'), symmetry=1, barrier=(31.7025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37968,'amu*angstrom^2'), symmetry=1, barrier=(31.7215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56402,0.0392126,4.183e-06,-5.20864e-08,2.89632e-11,10219.1,16.198], Tmin=(100,'K'), Tmax=(886.039,'K')), NASAPolynomial(coeffs=[19.5025,-0.00143147,4.69969e-06,-1.09226e-09,7.70088e-14,5456.87,-77.1099], Tmin=(886.039,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.1241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
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
    label = '[CH2][C]=CC=C=C(18817)',
    structure = SMILES('C=[C]C=C[C]=C'),
    E0 = (543.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63111,'amu*angstrom^2'), symmetry=1, barrier=(37.5025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62034,'amu*angstrom^2'), symmetry=1, barrier=(37.2548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759229,0.0600512,-5.80612e-05,2.91516e-08,-5.61568e-12,65516.4,20.7626], Tmin=(100,'K'), Tmax=(1428.69,'K')), NASAPolynomial(coeffs=[14.7311,0.0143872,-3.2453e-06,3.65966e-10,-1.74471e-14,62192.1,-49.2902], Tmin=(1428.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=C(O)C[C]=C(28583)',
    structure = SMILES('[CH2][C]=C(O)C[C]=C'),
    E0 = (465.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,369.125,369.126],'cm^-1')),
        HinderedRotor(inertia=(0.163139,'amu*angstrom^2'), symmetry=1, barrier=(15.7738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163141,'amu*angstrom^2'), symmetry=1, barrier=(15.7738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163139,'amu*angstrom^2'), symmetry=1, barrier=(15.7738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163139,'amu*angstrom^2'), symmetry=1, barrier=(15.7738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.485275,0.0695245,-6.19063e-05,2.38095e-08,-2.17398e-12,56101.9,28.5297], Tmin=(100,'K'), Tmax=(1000.22,'K')), NASAPolynomial(coeffs=[16.8903,0.0182222,-6.42003e-06,1.12381e-09,-7.72414e-14,52104.7,-54.1953], Tmin=(1000.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CC(O)=C[CH2](20296)',
    structure = SMILES('[CH2][C]=CC(O)=C[CH2]'),
    E0 = (281.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835691,0.0580437,-2.16502e-05,-2.25839e-08,1.59892e-11,34033,25.798], Tmin=(100,'K'), Tmax=(919.664,'K')), NASAPolynomial(coeffs=[17.1659,0.0172447,-4.40847e-06,6.57182e-10,-4.40493e-14,29751.1,-58.568], Tmin=(919.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC(O)[CH][C]=C(20300)',
    structure = SMILES('[CH]=CC(O)[CH][C]=C'),
    E0 = (445.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,388.895,389.455],'cm^-1')),
        HinderedRotor(inertia=(0.245568,'amu*angstrom^2'), symmetry=1, barrier=(26.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111417,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244999,'amu*angstrom^2'), symmetry=1, barrier=(26.2217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244898,'amu*angstrom^2'), symmetry=1, barrier=(26.2312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972877,0.0523075,-7.87135e-06,-3.2579e-08,1.73874e-11,53760.1,30.8486], Tmin=(100,'K'), Tmax=(980.407,'K')), NASAPolynomial(coeffs=[17.547,0.0180641,-6.5473e-06,1.24618e-09,-9.27981e-14,48906.1,-56.9652], Tmin=(980.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]C(O)[C]=C(28734)',
    structure = SMILES('[CH]=C[CH]C(O)[C]=C'),
    E0 = (445.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,388.895,389.455],'cm^-1')),
        HinderedRotor(inertia=(0.245568,'amu*angstrom^2'), symmetry=1, barrier=(26.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111417,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244999,'amu*angstrom^2'), symmetry=1, barrier=(26.2217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244898,'amu*angstrom^2'), symmetry=1, barrier=(26.2312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972877,0.0523075,-7.87135e-06,-3.2579e-08,1.73874e-11,53760.1,30.8486], Tmin=(100,'K'), Tmax=(980.407,'K')), NASAPolynomial(coeffs=[17.547,0.0180641,-6.5473e-06,1.24618e-09,-9.27981e-14,48906.1,-56.9652], Tmin=(980.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC([O])[C]=C(18290)',
    structure = SMILES('C=[C]CC([O])[C]=C'),
    E0 = (596.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,337.232,337.287,337.292,4000],'cm^-1')),
        HinderedRotor(inertia=(0.237396,'amu*angstrom^2'), symmetry=1, barrier=(19.1604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237331,'amu*angstrom^2'), symmetry=1, barrier=(19.1606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23734,'amu*angstrom^2'), symmetry=1, barrier=(19.1603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13791,0.0609972,-5.09129e-05,2.22672e-08,-3.99731e-12,71899.5,29.3365], Tmin=(100,'K'), Tmax=(1307.01,'K')), NASAPolynomial(coeffs=[12.368,0.0266285,-1.14694e-05,2.14824e-09,-1.49045e-13,68963.9,-27.849], Tmin=(1307.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C][CH]C([O])C=C(18289)',
    structure = SMILES('C=[C][CH]C([O])C=C'),
    E0 = (429.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,546.233,546.346,546.388,546.508],'cm^-1')),
        HinderedRotor(inertia=(0.12446,'amu*angstrom^2'), symmetry=1, barrier=(26.3541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12436,'amu*angstrom^2'), symmetry=1, barrier=(26.3462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124248,'amu*angstrom^2'), symmetry=1, barrier=(26.341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.5,'J/mol'), sigma=(6.3848,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.16 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20991,0.0447833,1.47304e-05,-5.50627e-08,2.49446e-11,51740.9,29.403], Tmin=(100,'K'), Tmax=(978.524,'K')), NASAPolynomial(coeffs=[17.2383,0.019111,-7e-06,1.35861e-09,-1.02782e-13,46696.3,-57.3247], Tmin=(978.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH2][C]=C(O)C=C[CH2](28824)',
    structure = SMILES('[CH2][C]=C(O)C=C[CH2]'),
    E0 = (281.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835691,0.0580437,-2.16502e-05,-2.25839e-08,1.59892e-11,34033,25.798], Tmin=(100,'K'), Tmax=(919.664,'K')), NASAPolynomial(coeffs=[17.1659,0.0172447,-4.40847e-06,6.57182e-10,-4.40493e-14,29751.1,-58.568], Tmin=(919.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CC=CCJ)"""),
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
    label = 'C=[C]C([O])[CH]C=C(20299)',
    structure = SMILES('C=[C]C([O])[CH]C=C'),
    E0 = (429.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,527.406,527.647,528.058,529.206],'cm^-1')),
        HinderedRotor(inertia=(0.127981,'amu*angstrom^2'), symmetry=1, barrier=(24.9895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126621,'amu*angstrom^2'), symmetry=1, barrier=(25.0487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.598911,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2099,0.0447834,1.47302e-05,-5.50624e-08,2.49444e-11,51740.9,29.403], Tmin=(100,'K'), Tmax=(978.525,'K')), NASAPolynomial(coeffs=[17.2383,0.019111,-6.99998e-06,1.35861e-09,-1.02782e-13,46696.3,-57.3248], Tmin=(978.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(O)C[C]=C(28584)',
    structure = SMILES('[CH]=[C]C(O)C[C]=C'),
    E0 = (613.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1670,1700,300,440,287.328,287.338],'cm^-1')),
        HinderedRotor(inertia=(0.0020418,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201096,'amu*angstrom^2'), symmetry=1, barrier=(11.7808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201068,'amu*angstrom^2'), symmetry=1, barrier=(11.7808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201066,'amu*angstrom^2'), symmetry=1, barrier=(11.7808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01939,0.0670476,-6.81197e-05,3.76162e-08,-8.52214e-12,73913.9,30.3635], Tmin=(100,'K'), Tmax=(1055.48,'K')), NASAPolynomial(coeffs=[11.3942,0.02773,-1.22434e-05,2.32358e-09,-1.62811e-13,71723.8,-20.2492], Tmin=(1055.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(O)[CH][CH][CH2](23689)',
    structure = SMILES('C#CC(O)[CH][CH][CH2]'),
    E0 = (528.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3000,3050,390,425,1340,1360,335,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2175,525,3615,1277.5,1000,180,964.726],'cm^-1')),
        HinderedRotor(inertia=(3.42526,'amu*angstrom^2'), symmetry=1, barrier=(78.7535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00398485,'amu*angstrom^2'), symmetry=1, barrier=(2.63186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0192643,'amu*angstrom^2'), symmetry=1, barrier=(12.7231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114461,'amu*angstrom^2'), symmetry=1, barrier=(2.63169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126762,'amu*angstrom^2'), symmetry=1, barrier=(78.7521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01047,0.063721,-6.26858e-05,3.43933e-08,-7.66383e-12,63692,33.1994], Tmin=(100,'K'), Tmax=(1084.31,'K')), NASAPolynomial(coeffs=[11.4564,0.0251862,-9.3779e-06,1.61799e-09,-1.07105e-13,61426.7,-18.0416], Tmin=(1084.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJC) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH]=[C][CH]C(O)C=C(20060)',
    structure = SMILES('[CH]=[C][CH]C(O)C=C'),
    E0 = (445.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,388.895,389.455],'cm^-1')),
        HinderedRotor(inertia=(0.245568,'amu*angstrom^2'), symmetry=1, barrier=(26.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111417,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244999,'amu*angstrom^2'), symmetry=1, barrier=(26.2217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244898,'amu*angstrom^2'), symmetry=1, barrier=(26.2312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972877,0.0523075,-7.87135e-06,-3.2579e-08,1.73874e-11,53760.1,30.8486], Tmin=(100,'K'), Tmax=(980.407,'K')), NASAPolynomial(coeffs=[17.547,0.0180641,-6.5473e-06,1.24618e-09,-9.27981e-14,48906.1,-56.9652], Tmin=(980.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC=[C][CH2](18813)',
    structure = SMILES('[CH2][C]=CC=[C][CH2]'),
    E0 = (733.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,708.396],'cm^-1')),
        HinderedRotor(inertia=(2.93563,'amu*angstrom^2'), symmetry=1, barrier=(67.4959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190931,'amu*angstrom^2'), symmetry=1, barrier=(67.4582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644865,'amu*angstrom^2'), symmetry=1, barrier=(67.3549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3375.74,'J/mol'), sigma=(5.79513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.28 K, Pc=39.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84928,0.0398838,-8.2793e-06,-1.72322e-08,9.68703e-12,88323.6,22.6963], Tmin=(100,'K'), Tmax=(960.212,'K')), NASAPolynomial(coeffs=[10.3621,0.0219425,-7.62276e-06,1.3151e-09,-8.95683e-14,85881.1,-22.2335], Tmin=(960.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C][CH]C([O])[C]=C(20303)',
    structure = SMILES('C=[C][CH]C([O])[C]=C'),
    E0 = (667.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,562.69,562.845,564.729,564.788],'cm^-1')),
        HinderedRotor(inertia=(0.112337,'amu*angstrom^2'), symmetry=1, barrier=(25.4303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113075,'amu*angstrom^2'), symmetry=1, barrier=(25.4158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113175,'amu*angstrom^2'), symmetry=1, barrier=(25.4255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17167,0.0494391,-1.02245e-05,-2.56226e-08,1.3928e-11,80344.4,30.0143], Tmin=(100,'K'), Tmax=(998.562,'K')), NASAPolynomial(coeffs=[16.1254,0.0184837,-7.20539e-06,1.39075e-09,-1.02723e-13,75914.9,-49.3332], Tmin=(998.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(667.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC(O)=[C][CH2](28825)',
    structure = SMILES('[CH2][C]=CC(O)=[C][CH2]'),
    E0 = (519.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,330.516],'cm^-1')),
        HinderedRotor(inertia=(0.988767,'amu*angstrom^2'), symmetry=1, barrier=(76.643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.988697,'amu*angstrom^2'), symmetry=1, barrier=(76.6425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287012,'amu*angstrom^2'), symmetry=1, barrier=(22.25,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.98872,'amu*angstrom^2'), symmetry=1, barrier=(76.6425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296285,0.0685847,-6.7107e-05,3.34062e-08,-6.37367e-12,62658.5,28.2092], Tmin=(100,'K'), Tmax=(1425.69,'K')), NASAPolynomial(coeffs=[17.1342,0.0147491,-3.52754e-06,4.31579e-10,-2.25599e-14,58527.5,-56.6455], Tmin=(1425.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=[C]C(O)[CH][C]=C(28826)',
    structure = SMILES('[CH]=[C]C(O)[CH][C]=C'),
    E0 = (683.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1670,1700,300,440,180.846,706.324],'cm^-1')),
        HinderedRotor(inertia=(0.0593021,'amu*angstrom^2'), symmetry=1, barrier=(21.0966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.49814,'amu*angstrom^2'), symmetry=1, barrier=(21.1264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.059499,'amu*angstrom^2'), symmetry=1, barrier=(21.1395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0579646,'amu*angstrom^2'), symmetry=1, barrier=(21.1196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.92829,0.0570369,-3.30747e-05,-2.83427e-09,6.25061e-12,82364,31.4828], Tmin=(100,'K'), Tmax=(1009.59,'K')), NASAPolynomial(coeffs=[16.4734,0.0173712,-6.71536e-06,1.26957e-09,-9.20193e-14,78107.8,-49.1953], Tmin=(1009.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(O)[C]=C(22023)',
    structure = SMILES('[CH]=[C][CH]C(O)[C]=C'),
    E0 = (683.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1670,1700,300,440,244.921,709.027],'cm^-1')),
        HinderedRotor(inertia=(0.0589693,'amu*angstrom^2'), symmetry=1, barrier=(21.1138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0589364,'amu*angstrom^2'), symmetry=1, barrier=(21.1166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0592513,'amu*angstrom^2'), symmetry=1, barrier=(21.1179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919193,'amu*angstrom^2'), symmetry=1, barrier=(21.1341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.92828,0.057037,-3.30752e-05,-2.83369e-09,6.25036e-12,82364,31.4828], Tmin=(100,'K'), Tmax=(1009.59,'K')), NASAPolynomial(coeffs=[16.4734,0.0173712,-6.71531e-06,1.26956e-09,-9.20185e-14,78107.8,-49.1956], Tmin=(1009.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(O)=C=C(28580)',
    structure = SMILES('C=[C]CC(O)=C=C'),
    E0 = (252.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.793493,'amu*angstrom^2'), symmetry=1, barrier=(18.244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794948,'amu*angstrom^2'), symmetry=1, barrier=(18.2774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793346,'amu*angstrom^2'), symmetry=1, barrier=(18.2406,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.567029,0.0669534,-5.43367e-05,1.58449e-08,6.59769e-13,30513.4,26.5357], Tmin=(100,'K'), Tmax=(989.107,'K')), NASAPolynomial(coeffs=[16.943,0.0179974,-6.28333e-06,1.10839e-09,-7.70806e-14,26429.1,-56.5599], Tmin=(989.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C(O)C=C=C(20306)',
    structure = SMILES('[CH2]C=C(O)C=C=C'),
    E0 = (102.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625392,0.0620265,-2.93286e-05,-1.71872e-08,1.44525e-11,12477.1,23.5744], Tmin=(100,'K'), Tmax=(929.219,'K')), NASAPolynomial(coeffs=[18.7889,0.0150132,-3.76158e-06,5.74945e-10,-4.0105e-14,7755.66,-69.9628], Tmin=(929.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
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
    label = '[CH2][C]=C[CH]O(19108)',
    structure = SMILES('[CH2][C]=C[CH]O'),
    E0 = (323.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3615,1277.5,1000,1685,370,3010,987.5,1337.5,450,1655,1023.06],'cm^-1')),
        HinderedRotor(inertia=(0.0695788,'amu*angstrom^2'), symmetry=1, barrier=(20.2457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.70952,'amu*angstrom^2'), symmetry=1, barrier=(108.281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0271971,'amu*angstrom^2'), symmetry=1, barrier=(20.2596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19882,0.0305185,2.49335e-06,-2.72236e-08,1.31482e-11,38992.9,20.5429], Tmin=(100,'K'), Tmax=(969.873,'K')), NASAPolynomial(coeffs=[11.6505,0.013335,-4.64238e-06,8.53957e-10,-6.22871e-14,36134.4,-30.0518], Tmin=(969.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C(O)[C]=C(28828)',
    structure = SMILES('[CH2][C]=[C]C(O)[C]=C'),
    E0 = (739.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1670,1685,1700,300,370,440,1380,1390,370,380,2900,435,321.604,1814.83],'cm^-1')),
        HinderedRotor(inertia=(0.173793,'amu*angstrom^2'), symmetry=1, barrier=(12.7556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173747,'amu*angstrom^2'), symmetry=1, barrier=(12.7558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28925,'amu*angstrom^2'), symmetry=1, barrier=(21.2349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173754,'amu*angstrom^2'), symmetry=1, barrier=(12.7557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11419,0.0679443,-8.21642e-05,5.66747e-08,-1.61952e-11,89091.2,29.7803], Tmin=(100,'K'), Tmax=(843.374,'K')), NASAPolynomial(coeffs=[9.21358,0.0295314,-1.38464e-05,2.67296e-09,-1.88069e-13,87725,-7.91518], Tmin=(843.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(739.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C#CC(O)[C]=C(28829)',
    structure = SMILES('[CH2]C#CC(O)[C]=C'),
    E0 = (428.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,2100,2250,500,550,1380,1390,370,380,2900,435,253.119,253.28],'cm^-1')),
        HinderedRotor(inertia=(0.00263106,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00263054,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519216,'amu*angstrom^2'), symmetry=1, barrier=(23.6278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13599,'amu*angstrom^2'), symmetry=1, barrier=(51.6889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21882,0.0589799,-5.46774e-05,2.74703e-08,-5.61989e-12,51690.9,28.4514], Tmin=(100,'K'), Tmax=(1170.66,'K')), NASAPolynomial(coeffs=[11.5573,0.0236546,-9.41385e-06,1.69363e-09,-1.15129e-13,49270.3,-23.0546], Tmin=(1170.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_S) + radical(Propargyl)"""),
)

species(
    label = '[CH2]C=[C]C(O)[C]=C(28830)',
    structure = SMILES('[CH2]C=[C]C(O)[C]=C'),
    E0 = (502.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.171096,'amu*angstrom^2'), symmetry=1, barrier=(3.93384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0399424,'amu*angstrom^2'), symmetry=1, barrier=(15.6744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.039899,'amu*angstrom^2'), symmetry=1, barrier=(15.681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681608,'amu*angstrom^2'), symmetry=1, barrier=(15.6715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1666,0.0632582,-5.76434e-05,2.83305e-08,-5.77743e-12,60486.7,29.108], Tmin=(100,'K'), Tmax=(1154.57,'K')), NASAPolynomial(coeffs=[11.1295,0.0287421,-1.28014e-05,2.43835e-09,-1.71066e-13,58186.1,-20.3895], Tmin=(1154.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C(O)C=C(20311)',
    structure = SMILES('[CH2][C]=[C]C(O)C=C'),
    E0 = (502.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.171096,'amu*angstrom^2'), symmetry=1, barrier=(3.93384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0399424,'amu*angstrom^2'), symmetry=1, barrier=(15.6744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.039899,'amu*angstrom^2'), symmetry=1, barrier=(15.681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681608,'amu*angstrom^2'), symmetry=1, barrier=(15.6715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1666,0.0632582,-5.76434e-05,2.83305e-08,-5.77743e-12,60486.7,29.108], Tmin=(100,'K'), Tmax=(1154.57,'K')), NASAPolynomial(coeffs=[11.1295,0.0287421,-1.28014e-05,2.43835e-09,-1.71066e-13,58186.1,-20.3895], Tmin=(1154.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)[C]=[C]C(28831)',
    structure = SMILES('C=[C]C(O)[C]=[C]C'),
    E0 = (588.416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1685,1700,300,370,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,274.508,274.514],'cm^-1')),
        HinderedRotor(inertia=(0.150699,'amu*angstrom^2'), symmetry=1, barrier=(8.05575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15061,'amu*angstrom^2'), symmetry=1, barrier=(8.05596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150492,'amu*angstrom^2'), symmetry=1, barrier=(8.05613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150658,'amu*angstrom^2'), symmetry=1, barrier=(8.05605,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.966288,0.0723884,-9.86249e-05,8.21706e-08,-2.8389e-11,70874,30.0182], Tmin=(100,'K'), Tmax=(778.833,'K')), NASAPolynomial(coeffs=[7.19338,0.0350948,-1.65684e-05,3.17448e-09,-2.20829e-13,70065.1,2.56688], Tmin=(778.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(588.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C(O)C=[C]C(28832)',
    structure = SMILES('[CH2][C]=C(O)C=[C]C'),
    E0 = (401.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475867,0.0721335,-7.19471e-05,3.59122e-08,-6.57781e-12,48447,26.0032], Tmin=(100,'K'), Tmax=(989.147,'K')), NASAPolynomial(coeffs=[15.4287,0.0204913,-7.01727e-06,1.17084e-09,-7.69486e-14,45057.1,-48.1548], Tmin=(989.147,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(Cds_S) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]C([O])C=[C]C(20314)',
    structure = SMILES('C=[C]C([O])C=[C]C'),
    E0 = (580.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,249.802,249.83,255.07],'cm^-1')),
        HinderedRotor(inertia=(0.00278138,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255041,'amu*angstrom^2'), symmetry=1, barrier=(10.95,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247292,'amu*angstrom^2'), symmetry=1, barrier=(10.9537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46378,0.0600631,-5.39107e-05,2.7904e-08,-6.301e-12,69958.1,27.578], Tmin=(100,'K'), Tmax=(1019.62,'K')), NASAPolynomial(coeffs=[8.0363,0.0342787,-1.59782e-05,3.10217e-09,-2.1981e-13,68617.8,-4.25831], Tmin=(1019.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(O)C=[C]C(28833)',
    structure = SMILES('[CH]=[C]C(O)C=[C]C'),
    E0 = (597.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,1670,1700,300,440,3010,987.5,1337.5,450,1655,234.329],'cm^-1')),
        HinderedRotor(inertia=(0.243285,'amu*angstrom^2'), symmetry=1, barrier=(9.46422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242891,'amu*angstrom^2'), symmetry=1, barrier=(9.46239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242368,'amu*angstrom^2'), symmetry=1, barrier=(9.46444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242646,'amu*angstrom^2'), symmetry=1, barrier=(9.45903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06431,0.0694765,-8.30372e-05,5.88695e-08,-1.75561e-11,71984.5,29.6082], Tmin=(100,'K'), Tmax=(805.277,'K')), NASAPolynomial(coeffs=[8.358,0.0332466,-1.55508e-05,2.99875e-09,-2.10772e-13,70809.8,-3.99995], Tmin=(805.277,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC(O)=C=C(28834)',
    structure = SMILES('[CH2]C=CC(O)=C=C'),
    E0 = (102.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625392,0.0620265,-2.93286e-05,-1.71872e-08,1.44525e-11,12477.1,23.5744], Tmin=(100,'K'), Tmax=(929.219,'K')), NASAPolynomial(coeffs=[18.7889,0.0150132,-3.76158e-06,5.74945e-10,-4.0105e-14,7755.66,-69.9628], Tmin=(929.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
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
    label = 'C=C1[CH]C(O)C1=C(28811)',
    structure = SMILES('C=C1[CH]C(O)C1=C'),
    E0 = (83.4951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62092,0.0267869,7.7368e-05,-1.2783e-07,5.28281e-11,10151.1,21.8807], Tmin=(100,'K'), Tmax=(947.253,'K')), NASAPolynomial(coeffs=[20.228,0.0128022,-2.76397e-06,5.47877e-10,-5.08336e-14,3728.31,-82.1747], Tmin=(947.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.4951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=CCJC(O)C=C)"""),
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
    E0 = (436.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (549.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (636.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (670.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (684.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (602.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (572.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (623.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (599.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (551.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (600.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (672.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (578.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (561.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (758.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (586.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (822.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (837.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (755.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (762.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (885.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (731.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (895.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (895.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (464.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (459.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (439.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (835.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (958.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (951.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (655.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (537.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (706.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (677.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (750.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (567.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (640.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (729.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (514.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (531.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (445.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['C=C=CO(12571)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', 'C=[C]C=C(O)[C]=C(28823)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C#CC(O)[CH][C]=C(24342)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]=C=CC(O)[C]=C(27788)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=C(584)', 'C=[C]C=CO(19110)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C=CO(12571)', '[CH][C]=C(18825)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.00817,'m^3/(mol*s)'), n=1.99965, Ea=(13.9682,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;YJ] for rate rule [Cds_Ca;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(D)(132)', '[CH2][C]=CC=C=C(18817)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.73723,'m^3/(mol*s)'), n=1.84633, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;OJ_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -7.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]=C(O)C[C]=C(28583)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['[CH2][C]=CC(O)=C[CH2](20296)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=CC(O)[CH][C]=C(20300)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C[CH]C(O)[C]=C(28734)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.3587e+07,'s^-1'), n=1.74167, Ea=(154.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]CC([O])[C]=C(18290)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['C=[C][CH]C([O])C=C(18289)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['[CH2][C]=C(O)C=C[CH2](28824)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.30234e+06,'s^-1'), n=1.68744, Ea=(125.264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]CC(O)[C]=C(20061)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/O)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['C=[C]C([O])[CH]C=C(20299)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.13341e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_Cd;XH_out] for rate rule [R4HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C(O)C[C]=C(28584)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.73633e+06,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#CC(O)[CH][CH][CH2](23689)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5Hall;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C][CH]C(O)C=C(20060)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5Hall;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['OH(D)(132)', '[CH2][C]=CC=[C][CH2](18813)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.10333e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', 'C=[C][CH]C([O])[C]=C(20303)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.00518e+06,'m^3/(mol*s)'), n=0.282325, Ea=(7.09479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH2][C]=CC(O)=[C][CH2](28825)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH]=[C]C(O)[CH][C]=C(28826)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH]=[C][CH]C(O)[C]=C(22023)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['C=[C]CC(O)=C=C(28580)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['[CH2]C=C(O)C=C=C(20306)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['C=[C]C1C(=C)C1O(28827)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C]=CO(18753)', '[CH][C]=C(18825)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[C]=C(584)', '[CH2][C]=C[CH]O(19108)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(8)', '[CH2][C]=[C]C(O)[C]=C(28828)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(8)', '[CH2]C#CC(O)[C]=C(28829)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]=CO(18753)', 'C#C[CH2](17441)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C=[C]C(O)[C]=C(28830)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][C]=[C]C(O)C=C(20311)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_Cs;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C]C(O)[C]=[C]C(28831)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['[CH2][C]=C(O)C=[C]C(28832)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C]C([O])C=[C]C(20314)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C]C(O)C=[C]C(28833)'],
    products = ['C=[C][CH]C(O)[C]=C(20297)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['[CH2]C=CC(O)=C=C(28834)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['[CH2]C([CH]O)=C[C]=C(20132)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C][CH]C(O)[C]=C(20297)'],
    products = ['C=C1[CH]C(O)C1=C(28811)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

network(
    label = '5230',
    isomers = [
        'C=[C][CH]C(O)[C]=C(20297)',
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
    label = '5230',
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

