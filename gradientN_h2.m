function varargout=gradientN_h2(I,type,pasxy)
[w,d] = fspecialSeperable(type);
padOpt = 'symmetric';
%w=1
if length(size(I))==2
    if length(w) > 1

        Ix = -ndnanfilter(I,w',[],[],[],{'symmetric'},1);
        Iy = -ndnanfilter(I,w,[],[],[],{'symmetric'},1);

    else
        Ix = -I; Iy = -I;
    end
 [Ix,Iy] = derivnanfilter2(Ix,Iy,pasxy,d);
% Ix=-Ix;
% Iy=-Iy;
 
else
    if length(w) > 1
        Ix = ndnanfilter(I,w',[],[],[],{'symmetric'},1); Ix = ndnanfilter(Ix,reshape(w,1,1,[]),[],[],[],{'symmetric'},1);
        Iy = ndnanfilter(I,w,[],[],[],{'symmetric'},1);  Iy = ndnanfilter(Iy,reshape(w,1,1,[]),[],[],[],{'symmetric'},1);
        Iz = ndnanfilter(I,w,[],[],[],{'symmetric'},1);  Iz = ndnanfilter(Iz,w',[],[],[],{'symmetric'},1);
    size(Iz)
    else
        Ix = -I; Iy = -I; Iz = -I;
    end
    
[Ix,Iy,Iz]=derivnanfilter3(Ix,Iy,Iz,pasxy,d);
Ix=Ix;
Iy=-Iy;
Iz=-Iz;
    varargout{3} = Iz;
end

varargout{1} = Ix; varargout{2} = Iy;
function [w,d] = fspecialSeperable(type)

sumFlag01 = 1;
switch lower(type)
    case 'fb', w = 1; d = [1 -1]; sumFlag01 = 0;
     case 'prewitt', w = [1 1 1]; d = [1 0 -1];
    case 'sobel', w = [1 2 1]; d = [1 0 -1];
    case 'scharr', w = [3 10 3]; d = [1 0 -1];
    case 'stencil', w = 1; d = [-1 8 0 -8 1]/12; sumFlag01 = 0;
    otherwise
        if strcmpi(type(1:7),'optimal')
            sumFlag01 = 0;
            switch str2double(type(8:end))
                case 5
                    w = [0.0376593171958139,0.249153396177344,0.426374573253684,0.249153396177344,0.0376593171958139];
                    d = [0.109603762960256,0.276690988455550,0,-0.276690988455550,-0.109603762960256];
             
                case 7
                    w = [0.00541196740974425,0.0695905825057286,0.244559723794791,0.360875452579473,0.244559723794791,0.0695905825057286,0.00541196740974425];
                    d = [0.0194786630434688,0.123914729925000,0.193554838845068,0,-0.193554838845068,-0.123914729925000,-0.0194786630434688];
                case 9
                    w = [0.000737598362084457,0.0155298478667911,0.0902598182960924,0.234469365350285,0.318006740249494,0.234469365350285,0.0902598182960924,0.0155298478667911,0.000737598362084457];
                    d = [0.00303163095459785,0.0352414678518254,0.118879484725614,0.144382960330377,0,-0.144382960330377,-0.118879484725614,-0.0352414678518254,-0.00303163095459785];
                case 11
                    w = [9.73660046309032e-05,0.00304151665483548,0.0261779198435670,0.103249150042758,0.223724696218449,0.287418702471518,0.223724696218449,0.103249150042758,0.0261779198435670,0.00304151665483548,9.73660046309032e-05];
                    d = [0.000439905393637534,0.00808328124581577,0.0449782033510644,0.108830913335368,0.112870645898035,0,-0.112870645898035,-0.108830913335368,-0.0449782033510644,-0.00808328124581577,-0.000439905393637534];
                case 13
                    w = [1.49541584976590e-05,0.000608044910058278,0.00691160005872185,0.0369033951658625,0.112173742521741,0.212488560211325,0.261799405947588,0.212488560211325,0.112173742521741,0.0369033951658625,0.00691160005872185,0.000608044910058278,1.49541584976590e-05];
                    d = [7.14573394720322e-05,0.00177958343198347,0.0138717416627181,0.0506584387298514,0.0970421434222207,0.0891267488305157,0,-0.0891267488305157,-0.0970421434222207,-0.0506584387298514,-0.0138717416627181,-0.00177958343198347,-7.14573394720322e-05];
                case 15
                    w = [5.98860795150199e-06,0.000191445809195946,0.00198513495426173,0.0115248319580300,0.0440836748027682,0.114823905477918,0.203925200191027,0.246919636397694,0.203925200191027,0.114823905477918,0.0440836748027682,0.0115248319580300,0.00198513495426173,0.000191445809195946,5.98860795150199e-06];
                    d = [2.68289243388421e-05,0.000526939396634328,0.00397615629561318,0.0177322281254572,0.0506284057542549,0.0879544206261874,0.0780467781003990,0,-0.0780467781003990,-0.0879544206261874,-0.0506284057542549,-0.0177322281254572,-0.00397615629561318,-0.000526939396634328,-2.68289243388421e-05];
                case 17
                    w = [7.69392199800271e-07,3.43462804426615e-05,0.000460563932733138,0.00328564897335160,0.0153609347354413,0.0504415385506543,0.117907645586470,0.196238466137000,0.232540172823407,0.196238466137000,0.117907645586470,0.0504415385506543,0.0153609347354413,0.00328564897335160,0.000460563932733138,3.43462804426615e-05,7.69392199800271e-07];
                    d = [3.73783393466673e-06,0.000104597678828706,0.00102618747708425,0.00569710286194269,0.0209031813657965,0.0513851162956773,0.0801001019942813,0.0666269955155201,0,-0.0666269955155201,-0.0801001019942813,-0.0513851162956773,-0.0209031813657965,-0.00569710286194269,-0.00102618747708425,-0.000104597678828706,-3.73783393466673e-06];
                case 19
                    w = [1.52643475621529e-07,7.90799631627647e-06,0.000122710209565798,0.000998814542578885,0.00531312830085921,0.0202775623584220,0.0572293928770952,0.120087013816114,0.187334615968446,0.217257402574254,0.187334615968446,0.120087013816114,0.0572293928770952,0.0202775623584220,0.00531312830085921,0.000998814542578885,0.000122710209565798,7.90799631627647e-06,1.52643475621529e-07];
                    d = [7.63266760902856e-07,2.52762743204803e-05,0.000290470367850947,0.00185854511355811,0.00795173988567152,0.0240580568530862,0.0508916097992498,0.0712071454435055,0.0555263098864427,0,-0.0555263098864427,-0.0712071454435055,-0.0508916097992498,-0.0240580568530862,-0.00795173988567152,-0.00185854511355811,-0.000290470367850947,-2.52762743204803e-05,-7.63266760902856e-07];
                otherwise
                    w = [0.229878817299031,0.540242365401939,0.229878817299031];
                    d = [0.425286806086887,0,-0.425286806086887];
                 
            end
        else
            w = 1; d = [1 0 -1];
        end

end
        w = w/sum(abs(w));

if sumFlag01
    d = d/sum(abs(d));
end