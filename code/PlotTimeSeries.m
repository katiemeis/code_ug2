

classdef PlotTimeSeries < handle

    properties

        probesinputfile = '';
        datainputfile = '';
        data = {};
        name2probe = containers.Map;
        probe2row = containers.Map;

    end

    methods

        function self = PlotTimeSeries(genelist, datafile)
        %constuctor

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

        function plot_original_data(self, genes)

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
                plot(x,y)
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
                plot(x1,y1)
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



