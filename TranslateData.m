% Modify this function as needed to allow translation from any other format
% to that expected by the HILO function. Mandatory fields are listed below
% as fields of the "d_out" structure

function d_out = TranslateData( d_in, varargin )
    ip = inputParser;
    ip.addRequired('d_in', @isstruct);
    ip.addParamValue('crop', 1, @isscalar);
    ip.parse(d_in, varargin{:});
    n = ip.Results.crop;
    
    d_out.fname     = d_in.fname;
    d_out.headers 	= d_in.headers;
    d_out.n         = d_in.n;
    d_out.m         = d_in.m;
    d_out.Vg        = d_in.Vgs(1:end-n);
    d_out.Chf       = d_in.Ch(1:end-n);
    d_out.Cqs       = d_in.Cq(1:end-n);
    if isfield(d_in, 'G')
        d_out.G     = d_in.G(1:end-n);
    else
        d_out.G     = zeros(size(d_out.Vgs));
    end
    if isfield(d_in, 'Qt')
        d_out.Qt    = d_in.Q_t(1:end-n);
    else
        d_out.Qt    = zeros(size(d_out.Vgs));
    end
end