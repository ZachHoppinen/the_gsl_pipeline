function s = isalpha_num(a)
%ISALPH_NUM returns True for alpha-numeric character, including '_'.
%
%   Copyright 1986-2010 The MathWorks, Inc.
%

error(nargchk(1, 1, nargin, 'struct'))
if ~ischar(a) || length(a)~=1
   ctrlMsgUtils.error('Ident:utility:isalphaNum1')
end

s = ~isempty(findstr(upper(a), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'));

% FILE END
