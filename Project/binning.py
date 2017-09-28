from numpy import *
import numpy as np
import random
from matplotlib import pyplot as plt

#data =loadtxt("HMC_Results.dat", unpack=True)
myarray = np.fromfile('HMC_Results.dat')

data = [1.413959,1.985235,1.744739,1.473122,1.436441,1.260630,1.664726,1.757111,2.208004,1.775983,1.908975,1.863935,1.875755,1.325478,2.122043,2.154155,1.905836,1.899651,1.630728,1.755186,1.728184,1.536091,1.868960,1.465597,2.267266,1.151611,1.498982,1.795573,1.774674,1.349825,1.809829,1.893947,1.693693,1.783495,1.840136,1.534353,1.661597,1.742603,1.943542,1.680768,1.600840,1.711752,1.667057,1.781991,1.670590,2.135384,1.834102,1.521280,1.564322,1.903888,1.910803,1.592002,1.594045,2.177130,1.666030,1.319627,2.163564,1.442664,2.094136,1.408079,1.497037,1.687624,1.567027,1.462262,1.529521,1.498426,1.292260,1.450073,1.615761,1.587050,1.839729,1.628071,1.624806,1.717996,2.086699,1.740979,1.616753,1.992338,2.124483,1.867802,1.268810,1.741784,1.365115,1.831881,2.048399,1.923777,1.915267,1.747559,1.439214,1.651138,1.753116,1.756469,1.483060,1.788670,1.556353,1.350765,1.570957,1.744908,1.483979,1.533237,1.613484,1.581169,1.267465,1.812995,1.735873,2.110312,1.783734,1.668570,1.721018,1.826188,1.730467,2.079891,1.795165,1.489390,2.118291,1.356519,1.658260,1.931170,1.940041,1.815141,1.591327,1.794306,1.779088,1.739942,1.889888,1.642940,1.599739,1.536854,1.583268,1.640562,1.719675,2.189097,1.804774,1.500554,1.964695,1.345518,1.772207,1.784098,1.488394,1.496643,1.515029,1.577802,1.788279,1.793205,1.756878,1.536251,1.626081,2.245190,1.697896,1.329121,1.896372,1.618826,1.426629,1.749589,1.773985,1.368582,1.956210,1.715463,1.551510,1.871153,1.247552,1.735752,1.714222,1.562342,1.358370,1.695042,1.554703,1.143316,1.987523,1.530240,1.736845,1.787819,1.736221,1.794645,1.900053,1.557279,1.727829,1.419579,1.583197,1.728996,1.853953,1.479975,2.068023,1.706213,1.912922,1.876619,1.856124,1.650357,1.614714,1.443135,1.923848,2.155073,1.535188,2.090775,1.669913,1.654748,1.864633,1.608879,1.629277,2.050249,1.765112,1.738354,1.715743,1.525924,1.310046,2.068186,2.510492,1.353192,1.666644,1.440436,1.837930,1.832356,1.651095,1.914921,2.020091,1.573773,1.413141,1.855562,1.613105,1.653936,1.473045,1.630039,1.823182,1.867241,1.988511,1.453911,1.569086,1.753832,1.428086,1.864511,1.488036,1.623487,1.696494,1.799270,1.829431,1.137930,1.435625,1.550382,1.644004,1.669522,1.590363,1.433560,1.931825,1.615990,1.589498,1.384402,1.609429,1.586588,2.007962,1.701426,1.600037,1.816618,1.849769,2.007857,1.776434,1.795386,1.798433,2.194510,1.868132,1.584005,1.613685,1.853157,2.055889,1.689042,1.384137,1.740451,1.724049,1.964480,1.402838,2.043165,1.553961,1.720591,1.295089,1.627082,1.939271,1.920089,1.678405,1.813298,1.679582,2.604576,1.700949,1.495642,1.769366,1.789577,1.580588,2.039802,1.661305,1.837668,1.573215,1.686984,1.625138,1.669111,2.055762,2.242105,1.996452,1.875845,1.754312,1.802536,2.196817,1.698463,1.548303,1.897148,1.412644,1.927922,1.479725,1.745621,1.394310,2.006859,1.184694,1.760123,1.741647,2.057468,1.981126,1.424603,1.535326,1.615105,1.725805,1.637839,1.485585,1.494685,1.456398,1.723192,1.857535,1.744161,1.600786,1.742021,1.287588,1.743462,1.582634,1.639209,1.744158,1.423565,1.776956,1.611030,1.620655,1.768608,1.526737,1.737376,1.618285,1.660758,1.631249,2.019246,1.488698,1.930184,1.645212,1.568940,2.063447,1.311465,1.735670,1.683652,1.737801,1.473875,1.560369,1.686154,1.361949,1.538325,1.721243,1.827100,2.393770,1.873632,1.506768,1.763470,1.414105,1.710935,1.938197,1.771946,2.046963,1.953771,1.670051,1.657481,1.667143,1.903692,1.504522,1.614523,1.464245,1.790294,1.570846,1.205554,1.378299,1.897909,1.761152,1.525603,1.770924,2.194107,1.932505,1.787206,1.604237,1.671087,1.545229,1.455610,1.542545,1.689495,1.104789,1.709348,1.733979,1.846883,1.691903,1.710651,1.574474,2.127614,1.898373,1.773366,1.792643,1.581236,1.449501,1.494514,1.588557,1.589142,1.832665,1.541816,2.368514,1.689734,1.813107,2.043133,1.901514,1.384897,1.267370,1.863214,1.364247,1.553428,1.753939,1.913528,1.439174,1.794926,1.998835,1.538542,1.684878,1.356156,1.953860,1.659303,1.688264,1.440963,1.851010,1.495117,1.513071,1.513642,1.694150,1.457481,1.695018,1.680301,1.631734,1.558987,1.737456,1.716057,1.594010,1.652344,1.681266,1.499000,1.453020,2.044897,1.504149,1.379606,1.784779,1.581002,1.604764,1.478344,1.822574,1.742793,1.337882,1.891503,1.546431,2.252790,1.826502,2.165555,1.758694,1.720632,1.727726,1.741553,1.593908,1.490992,1.520480,1.415486,1.788254,1.534885,1.866244,1.849094,1.400661,1.550918,1.729857,2.023969,2.018683,1.614327,1.642116,2.053481,1.141997,1.403954,1.752655,1.714627,1.756938,1.794926,1.624918,1.870088,1.805780,1.758688,1.481211,1.536381,1.883361,1.442776,1.407037,1.769750,1.748345,1.937199,1.739504,1.681033,1.918754,1.867463,2.153768,1.278549,1.775452,1.730132,1.627180,1.647286,1.548211,1.333707,1.640571,1.555619,1.651258,1.971168,1.642095,2.080669,1.466561,1.615527,2.035865,2.253670,1.740928,1.997864,1.554272,1.470323,1.309503,1.574798,1.443141,1.532138,1.615233,1.403603,1.670847,1.849971,1.725086,2.168987,1.516457,2.233298,1.872735,1.503353,1.415240,2.442978,1.932851,1.884707,1.550931,1.702674,1.625111,1.747608,1.776642,1.836767,1.592008,1.901271,1.737699,1.953118,1.659402,1.796163,1.423702,1.661079,1.445398,1.544984,1.928717,1.957199,2.278741,1.442069,1.438394,1.681727,1.667223,2.088279,1.955155,2.043724,1.761391,1.346550,1.552602,1.765384,1.585600,1.427464,1.936125,1.441538,1.818170,1.762560,1.603235,1.782729,1.286207,1.849165,2.221458,2.180401,1.178982,1.617130,1.716399,1.846197,2.204635,1.542762,1.509393,2.032652,2.227820,1.891494,1.557705,1.715200,1.758859,1.254037,1.876686,2.263313,2.253703,1.924555,1.308258,1.655429,1.720210,1.773904,1.462715,1.900013,1.813260,1.457572,1.610932,1.678934,1.503061,1.956616,1.483677,1.740921,1.624582,1.476034,1.791542,1.554042,1.402644,2.105089,1.803520,1.773050,1.562977,1.836117,1.392188,1.662706,1.564207,1.723373,1.882048,2.466505,1.852144,1.475305,1.877861,1.721840,1.642307,1.711301,1.793811,1.900357,1.688170,1.587932,2.081255,2.002847,1.719500,1.403817,1.462334,1.892849,1.510633,1.370205,1.906125,1.942154,1.663095,1.450588,1.595642,1.959988,2.010131,1.633022,1.946225,1.840359,1.946052,1.918156,1.685692,1.636955,1.838921,1.687374,1.603223,1.830994,1.456554,1.800066,1.726992,1.818255,1.427227,2.046625,1.640099,2.307244,2.090508,1.511097,1.837169,2.187496,1.963818,1.732971,1.891314,1.327527,2.719853,2.070519,1.905176,1.451483,1.694372,1.809739,1.796818,1.851014,1.550633,1.610015,1.684655,1.855138,1.836062,1.428042,1.916463,1.530567,1.664520,1.580236,1.326994,2.119508,1.810256,2.293223,1.483448,1.845992,1.845584,1.744247,1.590335,1.492969,1.493499,1.665389,1.677852,2.285606,2.012738,1.944976,1.429358,1.938691,1.594872,1.422128,1.411413,1.301383,1.752544,1.483838,1.682945,1.563484,1.637685,2.055248,1.750848,1.539030,1.550829,1.708421,1.401046,1.835459,1.692304,1.652100,1.799798,1.373058,1.533305,2.294237,1.435091,1.518406,1.825193,1.502844,1.741914,1.599210,1.564221,2.277142,1.655129,1.960608,1.632732,1.805830,1.426295,1.558978,1.861379,1.585185,1.804105,1.657410,1.650952,1.694666,1.653486,1.400470,1.687399,2.069913,1.324085,1.811101,1.401563,2.108013,1.698914,1.712951,1.425216,1.610817,1.456231,1.556310,1.893041,1.642837,1.645481,1.714304,1.442066,1.158867,1.792608,2.017737,1.705455,1.367023,2.034283,2.177258,1.303578,1.683897,1.857144,1.925996,1.746979,1.551897,1.593521,1.594387,1.457269,1.943181,1.996898,1.508806,1.912049,1.631943,1.767544,1.657720,1.322275,1.393792,2.023386,1.706481,1.834137,2.308180,1.663638,2.290234,1.751504,1.666743,1.851850,1.615252,2.123415,1.797129,1.759057,1.363916,1.851814,1.575154,2.065470,1.665457,1.764441,1.548183,1.866542,1.549911,1.568399,1.667241,1.462360,1.710822,1.717933,2.029925,1.824697,1.257584,1.431480,1.720008,1.518262,1.802824,2.364469,1.633371,1.464256,1.655724,1.449928,1.887198,1.606236,2.142794,1.528666,1.575637,1.262370,2.085022,1.244308,1.418250,1.844090,1.740634,1.818085,1.573326,2.098870,1.729110,2.050234,1.634598,1.493123,1.821832,1.725269,1.750994,2.409259,1.337269,1.816563,2.171357,1.765344,1.590721,1.907987,1.820147,1.681360,2.144158,1.710067,1.899619,1.825387,1.360944,1.665216,1.578759,1.766287,1.778538,1.282198,1.342074,1.943888,1.645287,1.405035,1.708542,1.587530,1.794689,1.620630,1.873110,1.789870,1.877293,1.602570,2.016866,1.629813,2.107819,1.647211,1.856502,1.495812,1.843291,1.574578,1.339299,1.958087,1.771810,1.694935,1.619334,1.650428,2.097163,2.199461,1.962136,1.776302,2.013787,2.014306,1.567444,1.365684,1.989643,1.944195,2.038470,1.651936,1.593473,1.766582,2.328418,1.205877,1.644391,1.880074,1.155893,1.280168,1.562071,1.614602,1.401164,1.797463,1.527413,1.601493,1.572123,1.448113,1.494942,1.513202,1.465795,1.282448,1.657370,1.462266,1.702034,1.728386,2.014167,1.631486,1.886218,1.633745,2.146295,1.282531,1.575130,1.609921,1.627755,1.711887,1.611564,1.993229,1.853495,2.063852,1.351569,1.630114,1.900859,1.678924,1.357241,1.862664,1.648197,1.790204,2.045073,1.659226,1.728472,1.847514,1.548007,1.957965,1.714194,2.135068,1.722676,1.672744,1.482380,1.518237,1.849719,1.650979,1.874666,1.399348,1.615185,1.485720,1.473191,1.646297,2.323979,1.980159,1.518595,1.857205,1.739175,1.538048,2.185730]

plt.hist(data, bins=50)
axes = plt.gca()
#axes.set_xlim([-0.5,0.5])
plt.show()

