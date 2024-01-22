function monthnumber = monthname2monthnumber(monthname)
% MONTHNAME2MONTHNUMBER converts a 3-character string
% 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
% 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', or 'DEC' to
% the corresponding integer numbers 1-12.
% Note: the function is case sensitive.
%
% Othmar Frey, 24.03.2004
%

switch monthname
	case 'JAN'
    	monthnumber = 1;
    case 'FEB'
    	monthnumber = 2; 
    case 'MAR'
    	monthnumber = 3;
    case 'APR'
    	monthnumber = 4;
    case 'MAY'
    	monthnumber = 5; 
    case 'JUN'
    	monthnumber = 6; 
    case 'JUL'
    	monthnumber = 7;
    case 'AUG'
    	monthnumber = 8;
    case 'SEP'
    	monthnumber = 9; 
    case 'OCT'
    	monthnumber = 10;
    case 'NOV'
    	monthnumber = 11;
    case 'DEC'
    	monthnumber = 12;
    otherwise
    	errortext = sprintf('monthname must be either JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, or DEC');
    	errordlg(errortext,'Error in function monthname2monthnumber','modal'); 
end
