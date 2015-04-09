       % ---------------------------------------------------------------------------------------------------
        function [exploredProp, exploredLen] = plotPathFollowingTimecourse_1back(this, threshDist, varargin)
            % function plotFollowingTimecourse(this)
            %
            % This function returns the proportion of the trail pixels that the mouse came within a given
            % distance of, for each video frame.  It's nice because it does show the bouts of following in a 
            % way that even penalizes the animal for stopping along the trail. It's a  clean way of measuring 
            % the completeness of his trail exploration.
            % Algorthmically, this is a similar problem to finding the following segments except reversed.  As
            % implemented it's SLOW.  Thought for speeding it up is to do it on the full video once, to get the
            % frames where the animal is within distance of ANY pixel, then do the incremental on only those frames.
            %
            % exploredLen returns a length in mm using mm_conv
            pb=1;
            if (sum(sum(this.exploredProp)) > 0)
                calculated = 1;
            else 
                calculated = 0;
            end
            if nargin > 2
                ah = varargin{1};
                if isempty(ah)
                    pb = 0;
                end
            else
                % plotting part
                figure;
                ah = axes; hold on;
            end
            if ~calculated
                nPaths = length(this.paths);
                exploredProp = zeros(this.nFrames, nPaths);
                exploredLen = zeros(this.nFrames, nPaths);
                this.makePathsSkel();
                for ii = 1:this.nFrames
                    if ~mod(ii, 100)
                        disp('Completed 100 frames');
                    end
                    np = this.nosePos(1:ii,:);
                    nn = ~isnan(np(:,1));
                    np = np(nn,:);
                    for trailNum = 1:nPaths
                        trailPos = this.paths(trailNum).PixelList;
                        
                        distm = ipdm(single(np), single(trailPos));
                        [trailDist, mini] = nanmin(distm, [], 1);
                        explored = find(trailDist <= threshDist);
                        npx = length(this.paths(trailNum).PixelIdxList);
                        exploredProp(ii,trailNum) = length(explored)/npx;
                        exploredLen(ii,trailNum) = length(explored);
                    end
                end
                this.makePathsFull();
                this.exploredProp = exploredProp;
                this.exploredLen = exploredLen * this.mm_conv;
            else
                exploredProp = this.exploredProp;
                exploredLen = this.exploredLen;
            end
            
            % plotting part
            if pb
                plot(this.times/1000, this.exploredLen(:,1), 'g','LineWidth', 2); hold on;
                plot(this.times/1000, this.exploredLen(:,2), 'r','LineWidth', 2); 
                xlabel('Time (sec)');
                ylabel('Length of the trail explored');
            end
            
        end