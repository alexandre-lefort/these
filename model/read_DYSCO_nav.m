function torpedo_ = read_DYSCO_nav( filename_ )
%
% READ A DYSCO .NAV MODEL DEFINITION FILE
% FILL A STRUCTURE IN RETURN
%
% 2010/10/12 @NC - SIREHNA

fid_ = fopen( filename_ , 'rt' );
while ~feof(fid_)
    line_=fgets(fid_);
    ibeg_=find( ~isspace(line_));
    if ~isempty(ibeg_) && line_(ibeg_(1)) ~= '#' && line_(ibeg_(1)) ~= '['
        ieq_=find( line_ == '=' );
        if ~isempty(ieq_)
            var_=line_(1:ieq_(1)-1);
            var_=var_( ~isspace(var_) );
            switch(var_)
                case {'1/J'}
                    var_='unSurJ';
                otherwise
            end
            
            try
                torpedo_.(var_)=eval(['[ ' , line_(ieq_(1)+1:end) , ' ];']);
                torpedo_.(upper(var_))=torpedo_.(var_);
                torpedo_.(lower(var_))=torpedo_.(var_);
            catch
            end
        end
    end
end
fclose(fid_);
