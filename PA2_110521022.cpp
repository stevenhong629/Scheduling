#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <vector>
#include <istream>
#include <sstream>
using namespace std;

struct Node {
	int number;
	string type;
	vector<struct Node*> to;
	vector<struct Node*> from;
};

struct Edge {
    int from;
    int to;
};

struct ScheduledNode {
    struct Node *node;
    int scheduledTime;
};

int cost(struct Node *node) {
	string type = (*node).type;
	if (type == "*") return 3;
	return 1;	// 剩下三種情況+, i, o都算1
}

bool isScheduled(vector<struct ScheduledNode> *schedule, struct Node *node) {  
	for (int i=0; i<schedule->size(); i++) {
		if ((*schedule)[i].node == node) return true;
	}
	return false;
}

int getScheduledTime(vector<struct ScheduledNode> *schedule, struct Node *node) {   
	for (int i=0; i<schedule->size(); i++) {
		if ((*schedule)[i].node == node) return (*schedule)[i].scheduledTime;
	}
	return -1;
}

void operationsOfTime(vector<struct ScheduledNode> *schedule, int time, vector<struct Node *> *nodes) {    //輸入一個schedule和time，算出這個schedule
	for (int i=0; i<schedule->size(); i++) {
		struct ScheduledNode *vi = &(*schedule)[i];
		int runningTime = cost(vi->node);
		int startTime = vi->scheduledTime;
		int endTime = startTime + runningTime;
		if ((startTime <= time) && (time < endTime)) {
			nodes->push_back(vi->node);
		}
	}
}

int countAdder(vector<struct ScheduledNode> *schedule, int latency) {
	int max = -1;
	for (int i=0; i<=latency; i++) {
		vector<struct Node *> nodes;
		operationsOfTime(schedule, i, &nodes);
		int count = 0;
		for (int j=0; j<nodes.size(); j++) {
			if (nodes[j]->type == "+") count++;
		}
		if (count > max) max = count;
	}
	return max;
}

int countMultiplier(vector<struct ScheduledNode> *schedule, int latency) {
	int max = -1;
	for (int i=0; i<=latency; i++) {
		vector<struct Node *> nodes;
		operationsOfTime(schedule, i, &nodes);
		int count = 0;
		for (int j=0; j<nodes.size(); j++) {
			if (nodes[j]->type == "*") count++;
		}
		if (count > max) max = count;
	}
	return max;
}


int totalTime(vector<struct ScheduledNode> *schedule) {
	int max = -1;
	for(int i=0; i<schedule->size(); i++){
	   	int ti = (*schedule)[i].scheduledTime;
		if (ti > max) max = ti;
	}
	return max;
	// 算出schedule裡面所有ScheduledNode的scheduledTime的max
}

void printScheduleInOutputFormat(vector<struct ScheduledNode> *schedule, int latency, string outname) {
    
	ofstream outfile(outname, ios::out);
    outfile << countAdder(schedule, latency) << endl;
    outfile << countMultiplier(schedule, latency) << endl;
	int sum = 0;
    for(int t=1; t<=totalTime(schedule); t++) {
        vector<struct Node *> currentNodes;
        operationsOfTime(schedule, t, &currentNodes);
        int count = 0;
        for (int i=0; i<currentNodes.size(); i++) {
        	if (currentNodes[i]->type != "o") {
        		if (count >= 1) outfile << " ";	
            	outfile << currentNodes[i]->number;
				
            	count++;
            }
        }
    	if (count > 0)
		{
			sum++;
			if (sum < latency)
			{	
				outfile << endl;
				
			}
		}
    }
}

void printSchedule(vector<struct ScheduledNode> *schedule) {
	for (int i=0; i<schedule->size(); i++) {
		struct ScheduledNode *scheduledNode = &(*schedule)[i];
		int time = scheduledNode->scheduledTime;
		int nodeNumber = scheduledNode->node->number;
		printf("t=%d: v%d\n", time, nodeNumber);
	}
}

void doSchedule(vector<struct ScheduledNode> *schedule, struct Node *node, int time) {
	(*schedule).push_back({node, time});
}

// 把nodes裡面所有node的number和type印出來
void printNodes(vector<struct Node>* nodes) {
	for(int i=0;i<(*nodes).size();i++){
		cout<<(*nodes)[i].number<<" "<<(*nodes)[i].type<<endl;
	}
}

void printEdges(vector<struct Edge>* edges) {
	for(int i=0; i<(*edges).size(); i++) {   
    	cout << (*edges)[i].from << " " << (*edges)[i].to << endl;
    }
}

// 用ASAP演算法，幫所有nodes做排程，存到asapSchedule裡面
void asap(vector<struct ScheduledNode> *asapSchedule, vector<struct Node> *nodes) {
	// 把type為i的node都排入schedule
	for (int i=0; i<nodes->size(); i++) {
		if ((*nodes)[i].type == "i") {
			struct Node *v0 = &(*nodes)[i];
			doSchedule(asapSchedule, v0, 0);							
		}
	}
   	while (asapSchedule->size() < nodes->size()) {						// 還有nodes沒排進schedule
		for (int i=0; i<nodes->size(); i++) {
			struct Node *vi = &(*nodes)[i];
			if (isScheduled(asapSchedule, vi) == false) {				// vi還沒被排進schedule的話
																		// 判斷vi前面的node是否都排進schedule了
				bool allVjAreScheduled = true;							// 先假設vi前面全部的vj都排進去了
				for (int j=0; j<vi->from.size(); j++) {
					struct Node *vj = vi->from[j];
					if (isScheduled(asapSchedule, vj) == false) {		// 只要有一個vj其實還沒排進去
						allVjAreScheduled = false;						// 就代表這個假設是錯的
						break;											// 就可以不用再檢查下去了
					}
				}
				// 如果都排進去了，才去算vi應該被排到的時間
				if (allVjAreScheduled == true) {
					int max = -1;
					for (int j=0; j<vi->from.size(); j++) {				// 看vi前面得所有vj
						struct Node *vj = vi->from[j];
						int tj = getScheduledTime(asapSchedule, vj);
						if (tj + cost(vj) > max) { 						// 如果這個vj的tj+dj > max
							max = tj + cost(vj);						// 就把這個tj+dj設為新的max
						}
					}														
					// 這時候max就是所有vi前面的vj中，dj+tj的最大值	
					int ti = max;										// 把max作為vi要排的時間
					doSchedule(asapSchedule, vi, ti);					// schedule vi at ti
				}

			}
		}
	}
}

// 用ALAP演算法，幫所有nodes做排程，存到asapSchedule裡面
void alap(vector<struct ScheduledNode> *alapSchedule, int latency ,vector<struct Node> *nodes) {
	// 把type為o的node都排入schedule

	for (int i=0; i<nodes->size(); i++) {
		if ((*nodes)[i].type == "o") {
			struct Node *vi = &(*nodes)[i];
			doSchedule(alapSchedule, vi, latency + 1);
		}
	}
    
   	while (alapSchedule->size() < nodes->size()) {						// 還有nodes沒排進schedule
		for (int i=0; i<nodes->size(); i++) {
			struct Node *vi = &(*nodes)[i];

			if (isScheduled(alapSchedule, vi) == false) {				// vi還沒被排進schedule的話
																		// 判斷vi後面的node是否都排進schedule了
				bool allVjAreScheduled = true;							// 先假設vi後面全部的vj都排進去了
				for (int j=0; j<vi->to.size(); j++) {
					struct Node *vj = vi->to[j];
					if (isScheduled(alapSchedule, vj) == false) {		// 只要有一個vj其實還沒排進去
						allVjAreScheduled = false;						// 就代表這個假設是錯的
						break;											// 就可以不用再檢查下去了
					}
				}
				// 如果都排進去了，才去算vi應該被排到的時間
				if (allVjAreScheduled == true) {
					int min = latency + 2;
					for (int j=0; j<vi->to.size(); j++) {				// 看vi後面得所有vj
						struct Node *vj = vi->to[j];
						int tj = getScheduledTime(alapSchedule, vj);
						if (tj - cost(vi) < min) { 						// 如果這個vj的tj-di > min
							min = tj - cost(vi);						// 就把這個tj-di設為新的,om
						}
					}														
					// 這時候max就是所有vi前面的vj中，dj+tj的最大值	
					int ti = min;										// 把min作為vi要排的時間
					doSchedule(alapSchedule, vi, ti);					// schedule vi at ti
				}

			}
		}
	}
}

struct NodeTimeFrame {
	struct Node* node;
	int left;
	int right;	// inclusive
};

int findScheduledTimeByNode(vector<struct ScheduledNode> *schedule, struct Node *node) {   // 
    for (int j=0; j<schedule->size(); j++) {
        if ((*schedule)[j].node == node) {
            return (*schedule)[j].scheduledTime;
        }
    }
    return -1;
}

void computeTimeFrame(vector<struct ScheduledNode> *asapSchedule,
					  vector<struct ScheduledNode> *alapSchedule,
					  vector<struct NodeTimeFrame> *timeFrames) {
    for (int i=0; i<asapSchedule->size(); i++) {
        struct NodeTimeFrame frame;
        frame.node = (*asapSchedule)[i].node;
        frame.left = (*asapSchedule)[i].scheduledTime;
        frame.right = findScheduledTimeByNode(alapSchedule, frame.node);
        timeFrames->push_back(frame);
    }
}

void timeFrameToSchedule(vector<struct NodeTimeFrame> *timeFrames,
                        vector<struct ScheduledNode> *schedule) {
    for (int i=0; i<timeFrames->size(); i++) {
        struct NodeTimeFrame frame = (*timeFrames)[i];
        struct ScheduledNode scheduledNode;
        scheduledNode.node = frame.node;
        scheduledNode.scheduledTime = frame.left;
        schedule->push_back(scheduledNode);
    }
}

bool schedulingIsFinished(vector<struct NodeTimeFrame> *timeFrames) {
	for (int i=0; i<timeFrames->size(); i++) {
		struct NodeTimeFrame* vi = &(*timeFrames)[i];
		if (vi->left != vi->right) return false;
	}
	return true;
}

void scheduleOperation(vector<struct NodeTimeFrame> *oldTimeFrames,
						struct Node* node, int time,
						vector<struct NodeTimeFrame> *newTimeFrames) {
	for (int i=0; i<oldTimeFrames->size(); i++) {
		struct NodeTimeFrame* oldVi = &(*oldTimeFrames)[i];
		struct NodeTimeFrame newVi;
		newVi.node = node;
		if (oldVi->node == node) {
			newVi.left = time;
			newVi.right = time;
		} else {
			newVi.left = oldVi->left;
			newVi.right = oldVi->right;
		}
		newTimeFrames->push_back(newVi);
	}
}

// double totalForce(vector<struct NodeTimeFrame> *timeFrames) {

// }

bool isGreaterThan(double a, double b) {
	int scale = 1000;
	int newA = (int) (a * scale);
	int newB = (int) (b * scale);
	return (newA > newB);
}

/* 
 *   1. 根據asap和alap的結果，找出self-force的公式
 *   2. 取所有可能性的self-force最小的
*/
// void forceDirectedScheduling(vector<struct ScheduledNode> *asapSchedule,
// 				vector<struct ScheduledNode> *alapSchedule,
// 				int latency,
// 				vector<struct ScheduledNode> *forceDirectedSchedule) {

// 	vector<struct NodeTimeFrame> timeFrames;
// 	computeTimeFrame(asapSchedule, alapSchedule, &timeFrames);
// 	while (schedulingIsFinished(&timeFrames) == true) {
// 		double maxGain = -4147483648;	 // 檢查一下double的下限是不是這個值
// 		struct Node *bestOperation;
// 		int bestStep;
// 		for (int i=0; i<timeFrames.size(); i++) {
// 			struct NodeTimeFrame *vi = &(timeFrames[i]);
// 			if (vi->left == vi->right) continue;
// 			for (int j=vi->left; j<=vi->right; j++) {
// 				vector<struct NodeTimeFrame> workingFrames;
// 				scheduleOperation(&timeFrames, vi->node, j, &workingFrames);
// 				// double gain = totalForce(&timeFrames) - totalForce(&workingFrames);
// 				if (isGreaterThan(gain, maxGain)) {
// 					maxGain = gain;
// 					bestOperation = vi->node;
// 					bestStep = j;
// 				}
// 			}
// 		}
// 		vector<struct NodeTimeFrame> newFrames;
// 		scheduleOperation(&timeFrames, bestOperation, bestStep, &newFrames);
// 		timeFrames = newFrames;
// 	}
//     timeFrameToSchedule(&timeFrames, forceDirectedSchedule);
// }



void readInput(int *latency, vector<struct Node> *nodes, vector<struct Edge> *edges, string filename) {

    ifstream intestcase(filename, ios::in);
    
    
    struct Node newNode;
    struct Edge newEdge;
    int numberofNodes;
    int numberofEdges;
    string ss;
    /*while(getline(intestcase,ss))*/
    for (int i=0; i<3; i++) {
    	getline(intestcase, ss);
	    if (i == 0) *latency = stoi(ss);
	    else if (i == 1) numberofNodes = stoi(ss);
	    else numberofEdges = stoi(ss);
    }

 
    for (int j=0;j<numberofNodes;j++) {	 
  		string temp;
   		getline(intestcase, ss);
   		stringstream atr(ss);
   		atr >> temp;
   		newNode.number = stoi(temp);
   		atr >> temp;
   		newNode.type = temp;
  		(*nodes).push_back(newNode);
    }
    
    
    for(int i=0; i<numberofEdges; i++) {
		string temp1;
		getline(intestcase, ss);
		stringstream btr(ss);
		btr >> temp1;
		newEdge.from = stoi(temp1);
		btr >> temp1;
		newEdge.to = stoi(temp1);
		(*edges).push_back(newEdge);
		struct Node *to = &(*nodes)[newEdge.to - 1];
		struct Node *from = &(*nodes)[newEdge.from - 1];
		to->from.push_back(from);
		from->to.push_back(to);
    }
}



int main(int argc, char *argv[]) {
    
	
    // 把input資料存到nodes和edges兩個vector
    int latency;
    string filename, outname;
    filename = argv[1];
	outname = argv[2];
    
    vector<struct Node> nodes;
    vector<struct Edge> edges;
    readInput(&latency, &nodes, &edges, filename);
    
    // printNodes(&nodes);
    // printEdges(&edges);
	
    // 用ASAP計算出nodes的asapSchedule
	vector<struct ScheduledNode> asapSchedule;
	asap(&asapSchedule, &nodes);
    printScheduleInOutputFormat(&asapSchedule, latency, outname);
	// printSchedule(&asapSchedule);

    // 用ALAP計算出nodes的alapSchedule
	// vector<struct ScheduledNode> alapSchedule;
	// alap(&alapSchedule, latency, &nodes);
	// printSchedule(&alapSchedule);

	// vector<struct ScheduledNode> forceDirectedSchedule;
	// forceDirectedScheduling(&asapSchedule, &alapSchedule, latency, &forceDirectedSchedule);
    // printScheduleInOutputFormat(&forceDirectedSchedule, latency);
	// printSchedule(&forceDirectedSchedule);
	// printcountAdder;
	
}