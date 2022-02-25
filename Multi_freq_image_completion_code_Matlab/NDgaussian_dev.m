function dv_str = NDgaussian_dev(D,fsize,sigma)

	%filter size = 2 * fsize + 1
	filter_range = single(-fsize:fsize);
    
	% Create 1D gaussian distribution
	fil_space = 1/sigma/sqrt(2*pi)*exp(-0.5*filter_range.^2/sigma^2);
	fil_space = single(mat2gray(fil_space'));
	fil_temp = fil_space;
   
    % Extend it to N-D
	for i = 2:D
		fil_temp2 = fil_space;
		fil_space = fil_space*fil_temp(1);
		for j = 2:length(fil_temp)
			fil_space = cat(i,fil_space,fil_temp2*fil_temp(j));
		end
		fil_space = single(mat2gray(fil_space));
	end

	dv_fil = single(zeros([size(fil_space),D]));
	dv_str = '[';
	space_str = [];
	for i = 1:D;
		space_str = strcat(space_str,':,');
	end

	for i = 1:D
		dv_str = strcat(dv_str,'dv_fil(',space_str,num2str(i),')');
		if i ~= D
			dv_str = strcat(dv_str,',');
		else
			dv_str = strcat(dv_str,']');
		end
	end

	eval(sprintf(strcat(dv_str,' = gradient(fil_space);')));
end

