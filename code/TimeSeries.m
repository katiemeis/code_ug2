

classdef TimeSeries < handle

    properties

        probesinputfile = '';
        datainputfile = '';
        data = {};
        name2probe = containers.Map;
        probe2row = containers.Map;
        C = {};
        M = [];

    end

    methods

        function self = TimeSeries(genelist, datafile)
        %constuctor
        %genelist is the complete set of genes and probes
        %datafile is the full original time series matrix (not residuals)

            %read lines into cell array and map name to probe ID
            self.probesinputfile = genelist;
            fid = fopen(genelist, 'r');
            probes = textscan(fid, '%s%s%s%s%s%s%s', 'Delimiter', '\t');

            key = {probes{1,3}{:,1}};
            value = {probes{1,1}{:,1}};
            self.name2probe = containers.Map(key, value); 

            %read time sreies data file
            self.datainputfile = datafile;
            fid2 = fopen(datafile, 'r');
            fs = ['%q', repmat('%f', 1, 177)];
            self.data = textscan(fid2, fs, 'Delimiter', ',');

            %map probe name to row number/index
            key = self.data{1,1};
            index = [1:length(self.data{1,1})];
            self.probe2row = containers.Map(key, index);

        end

	function make_tweak_section(self, hemgenelist, adjfile, edgeweights)
        %makes T matrix, formats the tweak section and adds it to the stub file
        %hemgenelist only contains genes/probes you are using
        %adjfile contains only connections you want to model

            error(nargchk(3, 4, nargin));
            if nargin == 3
                edgeweights = 0;
            end

	    fid3 = fopen(hemgenelist);
	    self.C = textscan(fid3, '%s%s%s', 'Delimiter', '\t');

            %map gene name to index
	    key = {self.C{3}{:}};
            value = [1:length(key)];
            hashtable = containers.Map(key, value);

            self.M = zeros(length(hashtable.keys));

            %read adjfile and add 1 to appropriate place in M if connection exists
            fid4 = fopen(adjfile);
            line = fgetl(fid4);

            while ischar(line)

	        linecell = textscan(line, '%s');

                if edgeweights
                    self.M(hashtable(linecell{1}{1}), hashtable(linecell{1}{2}))=linecell{1}{3};
                else
		    self.M(hashtable(linecell{1}{1}), hashtable(linecell{1}{2}))=1;
                    self.M(hashtable(linecell{1}{2}), hashtable(linecell{1}{1}))=1;
                end

                line = fgetl(fid4);

	    end

            ones_row = repmat(1, 1, length(self.M));
            zeros_row = repmat(0, 1,length(self.M));

            %print tweak section to stub in proper format
            fid5 = fopen('teststub.txt', 'a+');
            f_string = [repmat('%d ', 1, length(self.M)), '\n'];
            fprintf(fid5, '%s\n', '$tweak');
            fprintf(fid5, '%s\n', 'promoter_strengths:');
            fprintf(fid5, f_string, ones_row);
            fprintf(fid5, '%s\n', 'genetic_interconnect_matrix:');
            fprintf(fid5, f_string, self.M);
            fprintf(fid5, '%s\n', 'external_input_strength:');
            fprintf(fid5, '%s\n', 'maternal_connection_strengths:');
            fprintf(fid5, f_string, ones_row);
            fprintf(fid5, '%s\n', 'promoter_thresholds:');
            fprintf(fid5, f_string, ones_row);
            fprintf(fid5, '%s\n', 'diffusion_parameters:');
            fprintf(fid5, f_string, zeros_row);
            fprintf(fid5, '%s\n', 'protein_half_lives:');
            fprintf(fid5, f_string, ones_row);
            fprintf(fid5, '%s\n', 'translational_transcriptional_delays:');
            fprintf(fid5, f_string, zeros_row);
            fprintf(fid5, '%s\n\n', '$$');
            fclose(fid5);

	end

        function make_input_section(self)
        %replaces 1s in T matrix with 0.001s, formats input section and adds it to stub file

            self.M(self.M==1)=0.001;

            fid5 = fopen('teststub.txt', 'a+');
            f_string = [repmat('%d ', 1, length(self.M)), '\n'];
            f_string2 = [repmat('%3.3f ', 1, length(self.M)), '\n'];
            f_string3 = [repmat('%2.1f ', 1, length(self.M)), '\n'];
            fprintf(fid5, '%s\n', '$input');
            fprintf(fid5, '%s\n', 'promoter_strengths:');
            fprintf(fid5, f_string, repmat(15, 1, length(self.M)));
            fprintf(fid5, '%s\n', 'genetic_interconnect_matrix:');
            fprintf(fid5, f_string2, self.M);
            fprintf(fid5, '%s\n', 'external_input_strength:');
            fprintf(fid5, '%s\n', 'maternal_connection_strengths:');
            fprintf(fid5, f_string2, repmat(.001, 1, length(self.M)));
            fprintf(fid5, '%s\n', 'promoter_thresholds:');
            fprintf(fid5, f_string3, repmat(-2.5, 1, length(self.M)));
            fprintf(fid5, '%s\n', 'diffusion_parameters:');
            fprintf(fid5, f_string, repmat(0, 1, length(self.M)));
            fprintf(fid5, '%s\n', 'protein_half_lives:');
            fprintf(fid5, f_string, repmat(10, 1, length(self.M)));
            fprintf(fid5, '%s\n', 'translational_transcriptional_delays:');
            fprintf(fid5, f_string, repmat(0, 1, length(self.M)));
            fprintf(fid5, '%s\n\n', '$$');
            fclose(fid5);

        

        end

        function format_model_data(self)
        %formats original time series data for the genes we need into format for bias and facts section, adds them to the stub file

            probes = {self.C{1,1}{:,1}};
            numprobes = length(probes);
            exp_matrix = [];

            %pull out the data for the probes you're using and add them to a matrix
            for j=1:numprobes

                probe = probes{j};
                rownum = self.probe2row(probe);

                dmatrix = {self.data(1,2:end)};
                D = dmatrix{1,1};
                datamatrix = cell2mat(D);
                genevect = datamatrix(rownum, :);

                newgenevect = reshape(genevect, 3, []);
                avgs = sum(newgenevect, 1)./size(newgenevect, 1);

                exp_matrix = [exp_matrix;avgs];

            end

            %add time points and lin columns
            zero_time_data = exp_matrix(:, 1);
            exp_matrix_2 = [zero_time_data'; exp_matrix'];
            timepoints = [0 0:2:24 27:3:48 52:4:72 96 120 168 2:2:24 27:3:48 52:4:72 96 120 168];

            time_mat = [timepoints;exp_matrix_2'];
            time_mat_trans = time_mat';
            time_mat_ordered = sortrows(time_mat_trans, 1);
            time_mat_new = time_mat_ordered';
            lin = [8192, 8193];       %8192 is for eryth and 8193 is for neut
            lin_row = [repmat(lin, 1, 30)];
            hem_matrix = [lin_row;time_mat_new];
            modeldata = hem_matrix';

            %time points 0 go to bias and the rest to facts
            bias = modeldata(1:2, :);
            facts = modeldata(3:end, :);

            %make header to title the columns in stub file
            titles = {'lin', 't'};
            genes_array = self.C{:, 3};
            genes_header = genes_array';
            header = [titles genes_header];

            %print bias and facts to stub in proper format
            fid5 = fopen('teststub.txt', 'a+');
            format_string = ['%d %d ', repmat('%f ', 1, length(self.M)), '\n'];
            fprintf(fid5, '%s\n', '$bias_wt');
            format_string2 = [repmat('%s ', 1, length(header)), '\n'];
            fprintf(fid5, format_string2, header{:});
            fprintf(fid5, format_string, bias');
            fprintf(fid5, '%s\n\n', '$$');

            fprintf(fid5, '%s\n', '$facts_wt');
            fprintf(fid5, format_string2, header{:});
            fprintf(fid5, format_string, facts');
            fprintf(fid5, '%s\n\n', '$$');
            fclose(fid5);

        end


        function plot_original_data(self, genes)
        %genes is a cell array containing genes you want to plot {'blah','merp'}

            numgenes = length(genes);
            colorarray = hsv(numgenes);
            ph = zeros(1, numgenes); 
            ph1 = zeros(1, numgenes);

            for j=1:numgenes
                gene = genes{j};

                %find row number needed
                probeID = self.name2probe(gene);
                rownum = self.probe2row(probeID);

                dmatrix = {self.data(1,2:end)};
                C = dmatrix{1,1};
                datamatrix = cell2mat(C);
                plotvect = datamatrix(rownum, :);

                %reshape and average triplicates
                newplotvect = reshape(plotvect, 3, []);
                avgpoints = sum(newplotvect, 1)./size(newplotvect, 1);

                %plot erythroid conditions
                subplot(2,1,1);
                x = [0:2:24 27:3:48 52:4:72 96 120 168];
                y = avgpoints(1, 1:30);
                hold on
                ph(j) = plot(x,y);
                set(ph(j), 'Color', colorarray(j,:))
                xlabel('Time (hrs)');
                ylabel('Expression');
                titlestring1 = sprintf('Erythroid Conditions');
                title(titlestring1);

                %plot neutrophil conditions
                subplot(2,1,2);
                x1 = [2:2:24 27:3:48 52:4:72 96 120 168];
                y1 = avgpoints(1, 31:end);
                hold on
                ph1(j) = plot(x1,y1);
                set(ph1(j), 'Color', colorarray(j,:))
                xlabel('Time (hrs)');
                ylabel('Expression');
                titlestring2 = sprintf('Neutrophil Conditions');
                title(titlestring2);

            end %from loop	    

	    legend(ph, genes);
        end %from function

        function plot_residuals(self, gene1, gene2)
        %uses original dataset, calculates residuals within

            %find the row we want
            probeID1 = self.name2probe(gene1);
            rownum1 = self.probe2row(probeID1);

            probeID2 = self.name2probe(gene2);
            rownum2 = self.probe2row(probeID2);

            %extract data from rows
            dmatrix = {self.data(1,2:end)};
            C = dmatrix{1,1};
            datamatrix = cell2mat(C);
            plotvect1 = datamatrix(rownum1, :);
            plotvect2 = datamatrix(rownum2, :);

            %reshape data and calculate residuals
            newplotvect1 = reshape(plotvect1, 3, []);
            meanvals1 = mean(newplotvect1);
            residuals1 = newplotvect1 - repmat(meanvals1, 3, 1);
            resvec1 = residuals1(:)';

            newplotvect2 = reshape(plotvect2, 3, []);
            meanvals2 = mean(newplotvect2);
            residuals2 = newplotvect2 - repmat(meanvals2, 3, 1);
            resvec2 = residuals2(:)';
            corr(resvec1', resvec2')

            %make scatter plot of erythrocyte residuals
            subplot(2,1,1);
            eryres1 = resvec1(1, 1:90);
            eryres2 = resvec2(1, 1:90);
            scatter(eryres1, eryres2)
            corr(eryres1', eryres2')
            xlabel(gene1);
            ylabel(gene2);
            titlestring1 = sprintf('Erythrocyte Residuals for %s and %s', gene1, gene2);
            title(titlestring1)

            %make scatter plot of neutrophil residuals
            subplot(2,1,2);
            neutres1 = resvec1(1, 91:end);
            neutres2 = resvec2(1, 91:end);
            scatter(neutres1, neutres2)
            corr(neutres1', neutres2')
            xlabel(gene1);
            ylabel(gene2);
            titlestring2 = sprintf('Neutrophil Residuals for %s and %s', gene1, gene2);
            title(titlestring2)


        end %from function


    end %from methods

end %from class



