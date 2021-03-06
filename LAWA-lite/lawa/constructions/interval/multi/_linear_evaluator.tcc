/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008-2011 Sebastian Kestler, Kristina Steih,
                         Alexander Stippler, Mario Rometsch.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI__LINEAR_EVALUATOR_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI__LINEAR_EVALUATOR_TCC 1

namespace lawa {
    
template <typename T>
T
_linear_bspline_inner_evaluator0(T x, unsigned short deriv){
    
    T value = 0.0L;
    
    if(deriv == 0){
        if(0. <= x && x < 0.5){
            value = (3.4641016151377545870548926830117447338860L)*x;
        } else if(0.5 <= x && x < 1.){
            value = (-3.4641016151377545870548926830117447338860L)*x + 3.4641016151377545870548926830117447338860L;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(0. <= x && x < 0.5){
            value = 3.4641016151377545870548926830117447338860L;
        } else if(0.5 <= x && x < 1.){
            value = -3.4641016151377545870548926830117447338860L;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}

template <typename T>
T
_linear_bspline_inner_evaluator1(T x, unsigned short deriv){
    
    T value = 0.0L;
    
    if(deriv == 0){
        if(0. <= x && x < 0.25){
            value = (9.4629797113899456762541985381637921932640L)*x;
        } else if(0.25 <= x && x < 0.5){
            value = (-14.1339733763590835017512525300328614367100L)*x + 5.8992382719372572945013627670491634074930L;
        } else if(0.5 <= x && x < 0.75){
            value = (2.9930033951944218584046121068203924559160L)*x - 2.6642501138394953855765695513774635388190L;
        } else if(0.75 <= x && x < 1.){
            value = (1.6779902697747159670924418850486767875270L)*x - 1.6779902697747159670924418850486767875270L;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(0. <= x && x < 0.25){
            value = 9.4629797113899456762541985381637921932640L;
        } else if(0.25 <= x && x < 0.5){
            value = -14.1339733763590835017512525300328614367100L;
        } else if(0.5 <= x && x < 0.75){
            value = 2.9930033951944218584046121068203924559160L;
        } else if(0.75 <= x && x < 1.){
            value = 1.6779902697747159670924418850486767875270L;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}

template <typename T>
T
_linear_bspline_inner_evaluator2(T x, unsigned short deriv){
    
    T value = 0.0L;
    
    if(deriv == 0){
        if(-1. <= x && x < -0.75){
            value = (0.4916304516561805397389447527790771058323L)*x + 0.4916304516561805397389447527790771058323L;
        } else if(-0.75 <= x && x < -0.5){
            value = (-3.4414131615932637781726132694535397408260L)*x - 2.4581522582809026986947237638953855291620L;
        } else if(-0.5 <= x && x < -0.25){
            value = (5.9269630459312163216395289863490475362870L)*x + 2.2260358454813373512113473640059081093950L;
        } else if(-0.25 <= x && x < 0.){
            value = (5.7077820374748175634619933623480694058970L)*x + 2.1712405933672376616669634580056635767970L;
        } else if(0. <= x && x < 0.25){
            value = (-11.0041910145033147414897859810894623670300L)*x + 2.1712405933672376616669634580056635767970L;
        } else if(0.25 <= x && x < 0.5){
            value = (2.9003852491162789750734587115532740655990L)*x - 1.3049034725376607674738477151550205313600L;
        } else if(0.5 <= x && x < 0.75){
            value = (-0.6780160427622340269601143229008770067205L)*x + 0.4842971734015957335429388020720550048003L;
        } else if(0.75 <= x && x < 1.){
            value = (0.0968594346803191467085877604144110009601L)*x - 0.0968594346803191467085877604144110009601L;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(-1. <= x && x < -0.75){
            value = 0.4916304516561805397389447527790771058323L;
        } else if(-0.75 <= x && x < -0.5){
            value = -3.4414131615932637781726132694535397408260L;
        } else if(-0.5 <= x && x < -0.25){
            value = 5.9269630459312163216395289863490475362870L;
        } else if(-0.25 <= x && x < 0.){
            value = 5.7077820374748175634619933623480694058970L;
        } else if(0. <= x && x < 0.25){
            value = -11.0041910145033147414897859810894623670300L;
        } else if(0.25 <= x && x < 0.5){
            value = 2.9003852491162789750734587115532740655990L;
        } else if(0.5 <= x && x < 0.75){
            value = -0.6780160427622340269601143229008770067205L;
        } else if(0.75 <= x && x < 1.){
            value = 0.0968594346803191467085877604144110009601L;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}


template <typename T>
T
_linear_wavelet_inner_evaluator0(T x, unsigned short deriv){
    
    T value = 0.0L;
    
    if(deriv == 0){
        if(0. <= x && x < 0.125){
            value = (0.5634759370212272862047463938812607934943L)*x;
        } else if(0.125 <= x && x < 0.25){
            value = (-5.1744137027564234980642342748537193529190L)*x + 0.7172362049722063480336225835918725183016L;
        } else if(0.25 <= x && x < 0.375){
            value = (9.7231265391058825214648985956349589346080L)*x - 3.0071488554933701568486606340302970535800L;
        } else if(0.375 <= x && x < 0.5){
            value = (9.4033649088454592233905444791693058039090L)*x - 2.8872382441457114200707778403556771295680L;
        } else if(0.5 <= x && x < 0.625){
            value = (-41.0288792947516617950954024122206874593500L)*x + 22.3288838576528490891721956053393195020600L;
        } else if(0.625 <= x && x < 0.75){
            value = (35.2183939040239842909577942565959518389100L)*x - 25.3256618915819297146110523126710800593500L;
        } else if(0.75 <= x && x < 0.875){
            value = (-6.4770915973142827125638712794283686254010L)*x + 5.9459522344217705380301968393471602888830L;
        } else if(0.875 <= x && x < 1.){
            value = (-2.2279766941741853162944757587787019332540L)*x + 2.2279766941741853162944757587787019332540L;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(0. <= x && x < 0.125){
            value = 0.5634759370212272862047463938812607934943L;
        } else if(0.125 <= x && x < 0.25){
            value = -5.1744137027564234980642342748537193529190L;
        } else if(0.25 <= x && x < 0.375){
            value = 9.7231265391058825214648985956349589346080L;
        } else if(0.375 <= x && x < 0.5){
            value = 9.4033649088454592233905444791693058039090L;
        } else if(0.5 <= x && x < 0.625){
            value = -41.0288792947516617950954024122206874593500L;
        } else if(0.625 <= x && x < 0.75){
            value = 35.2183939040239842909577942565959518389100L;
        } else if(0.75 <= x && x < 0.875){
            value = -6.4770915973142827125638712794283686254010L;
        } else if(0.875 <= x && x < 1.){
            value = -2.2279766941741853162944757587787019332540L;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_inner_evaluator1(T x, unsigned short deriv){
    
    T value = 0.0L;
    
    if(deriv == 0){
        if(-1. <= x && x < -0.75){
            value = (-0.4993361296332690191591578950775659138297L)*x - 0.4993361296332690191591578950775659138297L;
        } else if(-0.75 <= x && x < -0.5){
            value = (3.4953529074328831341141052655429613968080L)*x + 2.4966806481663450957957894753878295691480L;
        } else if(-0.5 <= x && x < -0.375){
            value = (-5.0211883451055647232098570219394899807000L)*x - 1.7615899781028788328661916683533961196050L;
        } else if(-0.375 <= x && x < -0.25){
            value = (-13.0105664192378690297563833431805446019700L)*x - 4.7576067559024929478211390388187916025830L;
        } else if(-0.25 <= x && x < -0.125){
            value = (6.2424769945550990144743973814052685475870L)*x + 0.0556540975457490632365561423276616848068L;
        } else if(-0.125 <= x && x < 0.){
            value = (5.7972442141891065085819482427839750691320L)*x;
        } else if(0. <= x && x < 0.125){
            value = (21.6687520790632535178928661292080107659000L)*x;
        } else if(0.125 <= x && x < 0.25){
            value = (-33.0912598879590913999048339023412547009400L)*x + 6.8450014958777931147247125039436581833550L;
        } else if(0.25 <= x && x < 0.375){
            value = (8.3814659203751117219235494602605469893280L)*x - 3.5231799562057576657323833367067922392110L;
        } else if(0.375 <= x && x < 0.5){
            value = (5.3297950450297485437320459474674898215430L)*x - 2.3788033779512464739105695194093958012920L;
        } else if(0.5 <= x && x < 0.625){
            value = (-1.3351060079635963904587827868469625109060L)*x + 0.9536471485454259931848448477478303649327L;
        } else if(0.625 <= x && x < 0.75){
            value = (-1.3351060079635963904587827868469625109060L)*x + 0.9536471485454259931848448477478303649327L;
        } else if(0.75 <= x && x < 0.875){
            value = (0.1907294297090851986369689695495660729865L)*x - 0.1907294297090851986369689695495660729865L;
        } else if(0.875 <= x && x < 1.){
            value = (0.1907294297090851986369689695495660729865L)*x - 0.1907294297090851986369689695495660729865L;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(-1. <= x && x < -0.75){
            value = -0.4993361296332690191591578950775659138297L;
        } else if(-0.75 <= x && x < -0.5){
            value = 3.4953529074328831341141052655429613968080L;
        } else if(-0.5 <= x && x < -0.375){
            value = -5.0211883451055647232098570219394899807000L;
        } else if(-0.375 <= x && x < -0.25){
            value = -13.0105664192378690297563833431805446019700L;
        } else if(-0.25 <= x && x < -0.125){
            value = 6.2424769945550990144743973814052685475870L;
        } else if(-0.125 <= x && x < 0.){
            value = 5.7972442141891065085819482427839750691320L;
        } else if(0. <= x && x < 0.125){
            value = 21.6687520790632535178928661292080107659000L;
        } else if(0.125 <= x && x < 0.25){
            value = -33.0912598879590913999048339023412547009400L;
        } else if(0.25 <= x && x < 0.375){
            value = 8.3814659203751117219235494602605469893280L;
        } else if(0.375 <= x && x < 0.5){
            value = 5.3297950450297485437320459474674898215430L;
        } else if(0.5 <= x && x < 0.625){
            value = -1.3351060079635963904587827868469625109060L;
        } else if(0.625 <= x && x < 0.75){
            value = -1.3351060079635963904587827868469625109060L;
        } else if(0.75 <= x && x < 0.875){
            value = 0.1907294297090851986369689695495660729865L;
        } else if(0.875 <= x && x < 1.){
            value = 0.1907294297090851986369689695495660729865L;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_inner_evaluator2(T x, unsigned short deriv){
    
    T value = 0.0L;
    
    if(deriv == 0){
        if(-1. <= x && x < -0.75){
            value = (-0.4916304516561805397389447527790771058323L)*x - 0.4916304516561805397389447527790771058323L;
        } else if(-0.75 <= x && x < -0.5){
            value = (3.4414131615932637781726132694535397408260L)*x + 2.4581522582809026986947237638953855291620L;
        } else if(-0.5 <= x && x < -0.375){
            value = (-3.9604412393064941626837499752327391129570L)*x - 1.2427749421689762717334578584477538977300L;
        } else if(-0.375 <= x && x < -0.25){
            value = (-19.6926156923042714343299820641632064995900L)*x - 7.1423403620431427486007948917966791677180L;
        } else if(-0.25 <= x && x < -0.125){
            value = (18.0000701462500477230961225830481207392500L)*x + 2.2808310975954370407557312700061526419920L;
        } else if(-0.125 <= x && x < 0.){
            value = (17.1233461124244526903859800870442082176900L)*x + 2.1712405933672376616669634580056635767970L;
        } else if(0. <= x && x < 0.125){
            value = (-33.0125730435099442244693579432683871010900L)*x + 2.1712405933672376616669634580056635767970L;
        } else if(0.125 <= x && x < 0.25){
            value = (22.6057320109684306417836208273025586294300L)*x - 4.7810475384425591966146588883157046395170L;
        } else if(0.25 <= x && x < 0.375){
            value = (-5.6124494201652150829139160031567820924810L)*x + 2.2734978193408522345597253192991305409600L;
        } else if(0.375 <= x && x < 0.5){
            value = (-2.5129475103950023882391076698956300617590L)*x + 1.1111846031770224740566721943261985294390L;
        } else if(0.5 <= x && x < 0.75){
            value = (0.6780160427622340269601143229008770067205L)*x - 0.4842971734015957335429388020720550048003L;
        } else if(0.75 <= x && x < 1.){
            value = (-0.0968594346803191467085877604144110009601L)*x + 0.0968594346803191467085877604144110009601L;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(-1. <= x && x < -0.75){
            value = -0.4916304516561805397389447527790771058323L;
        } else if(-0.75 <= x && x < -0.5){
            value = 3.4414131615932637781726132694535397408260L;
        } else if(-0.5 <= x && x < -0.375){
            value = -3.9604412393064941626837499752327391129570L;
        } else if(-0.375 <= x && x < -0.25){
            value = -19.6926156923042714343299820641632064995900L;
        } else if(-0.25 <= x && x < -0.125){
            value = 18.0000701462500477230961225830481207392500L;
        } else if(-0.125 <= x && x < 0.){
            value = 17.1233461124244526903859800870442082176900L;
        } else if(0. <= x && x < 0.125){
            value = -33.0125730435099442244693579432683871010900L;
        } else if(0.125 <= x && x < 0.25){
            value = 22.6057320109684306417836208273025586294300L;
        } else if(0.25 <= x && x < 0.375){
            value = -5.6124494201652150829139160031567820924810L;
        } else if(0.375 <= x && x < 0.5){
            value = -2.5129475103950023882391076698956300617590L;
        } else if(0.5 <= x && x < 0.75){
            value = 0.6780160427622340269601143229008770067205L;
        } else if(0.75 <= x && x < 1.){
            value = -0.0968594346803191467085877604144110009601L;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}


template <typename T>
T
_linear_wavelet_left_evaluator0(T x, unsigned short deriv){
    
    T value = 0.0L;
    
    if(deriv == 0){
        if(0. <= x && x < 0.125){
            value = (-26.6780669545463704769957840721302651348500L)*x;
        } else if(0.125 <= x && x < 0.25){
            value = (40.7411946788691781931779339150611206681700L)*x - 8.4274077041769435837717147483989232253780L;
        } else if(0.25 <= x && x < 0.375){
            value = (-10.3190672072464912686638287266090549462900L)*x + 4.3376577673519737816887259120186206782380L;
        } else if(0.375 <= x && x < 0.5){
            value = (-6.5619205271492485137235386635877667409980L)*x + 2.9287277623155077485861171383856376012530L;
        } else if(0.5 <= x && x < 0.625){
            value = (1.6437516725425437052863769025718135898170L)*x - 1.1741083375303883609188406446941525641550L;
        } else if(0.625 <= x && x < 0.75){
            value = (1.6437516725425437052863769025718135898170L)*x - 1.1741083375303883609188406446941525641550L;
        } else if(0.75 <= x && x < 0.875){
            value = (-0.2348216675060776721837681289388305128309L)*x + 0.2348216675060776721837681289388305128309L;
        } else if(0.875 <= x && x < 1.){
            value = (-0.2348216675060776721837681289388305128309L)*x + 0.2348216675060776721837681289388305128309L;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(0. <= x && x < 0.125){
            value = -26.6780669545463704769957840721302651348500L;
        } else if(0.125 <= x && x < 0.25){
            value = 40.7411946788691781931779339150611206681700L;
        } else if(0.25 <= x && x < 0.375){
            value = -10.3190672072464912686638287266090549462900L;
        } else if(0.375 <= x && x < 0.5){
            value = -6.5619205271492485137235386635877667409980L;
        } else if(0.5 <= x && x < 0.625){
            value = 1.6437516725425437052863769025718135898170L;
        } else if(0.625 <= x && x < 0.75){
            value = 1.6437516725425437052863769025718135898170L;
        } else if(0.75 <= x && x < 0.875){
            value = -0.2348216675060776721837681289388305128309L;
        } else if(0.875 <= x && x < 1.){
            value = -0.2348216675060776721837681289388305128309L;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_right_evaluator0(T x, unsigned short deriv){
    
    T value = 0.0L;
    
    if(deriv == 0){
        if(0. <= x && x < 0.25){
            value = (0.8560009184273419833656060881065680156050L)*x;
        } else if(0.25 <= x && x < 0.5){
            value = (-5.9920064289913938835592426167459761092350L)*x + 1.7120018368546839667312121762131360312100L;
        } else if(0.5 <= x && x < 0.625){
            value = (8.6077124804963414864334215100410800297640L)*x - 5.5878576178891837182651198871803920382900L;
        } else if(0.625 <= x && x < 0.75){
            value = (22.3037271753338132202831189197461682794400L)*x - 14.1478668021626035519211807682460721943400L;
        } else if(0.75 <= x && x < 0.875){
            value = (-10.7013406871343986753835187180574150506900L)*x + 10.6059340946885553698287974601066153032600L;
        } else if(0.875 <= x && x < 1.){
            value = (-9.9380879475676522309457486544510170712590L)*x + 9.9380879475676522309457486544510170712590L;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(0. <= x && x < 0.25){
            value = 0.8560009184273419833656060881065680156050L;
        } else if(0.25 <= x && x < 0.5){
            value = -5.9920064289913938835592426167459761092350L;
        } else if(0.5 <= x && x < 0.625){
            value = 8.6077124804963414864334215100410800297640L;
        } else if(0.625 <= x && x < 0.75){
            value = 22.3037271753338132202831189197461682794400L;
        } else if(0.75 <= x && x < 0.875){
            value = -10.7013406871343986753835187180574150506900L;
        } else if(0.875 <= x && x < 1.){
            value = -9.9380879475676522309457486544510170712590L;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}


template <typename T>
T
_linear_refinement_inner_evaluator0(T x, unsigned short deriv){

    T value = 0.0L;

    if(deriv == 0){
        if(0. <= x && x < 1.){
            value = x;
        } else if(1. <= x && x < 2.){
            value = 2.-x;
        } else {
            value = 0.0L;
        }
    } else if (deriv == 1){
        if(0. <= x && x < 1.){
            value = 1.;
        } else if(1. <= x && x < 2.){
            value = -1.;
        } else {
            value = 0.0L;
        }
    } else {
        value = 0.0L;
    }
    return value;
}
} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI__LINEAR_EVALUATOR_TCC
