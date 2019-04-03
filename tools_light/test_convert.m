model = load_model('iver2.nav');
sym_model = symtbx_convert_struct_to_sym(model);

f= fields(sym_model);
n = length(f);

for ii=1:n
    assert(abs(model.(f{ii}) - eval(sym_model.(f{ii}))) < 1e-10);
end