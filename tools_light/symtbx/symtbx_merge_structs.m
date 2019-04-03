function res = symtbx_merge_structs(varargin)
    res = [];
    for ii=1:nargin
       s = varargin{ii};
       f = fields(s);
       for jj=1:length(f)
           if ((ii + jj > 2) && isfield(res,f{jj}))
              warning(['element ' f{jj} 'is ereased']);
           end
           res.(f{jj}) = s.(f{jj});
       end
    end
