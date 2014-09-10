#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <stdlib.h>

typedef struct mnode{
  int mval; //min/max value
  //  int l; //left boundary
  //  int r; //right boundary
  struct mnode *prev;
  struct mnode *next;
  struct ynode *ly;
  //  bool released;
} mnode_t;


typedef struct ynode{
  int yval; //y value
  mnode_t *uptr;//pointer to its unode
  mnode_t *lptr;//pointer to its lnode
  struct ynode *prev;
  struct ynode *next;
} ynode_t;


mnode_t *ulist_head, *ulist_tail;
mnode_t *llist_head, *llist_tail;
ynode_t *ylist_head, *ylist_tail;

#define MAXMNODES 10000
#define MAXYNODES 10000

mnode_t mnodepool[MAXMNODES];
ynode_t ynodepool[MAXYNODES];

int mpoolptr;
int ypoolptr;


#define RECURSIVE 1


using namespace std;

vector <vector<int> > alignment;
vector <int> fcounts, ecounts;
vector <int> u, l; //upper/lower bounding the alignments column-wise
vector <int> lparens, rparens;

int n, m;


mnode_t *getmnode(){
  return &(mnodepool[mpoolptr++]);
}

ynode_t *getynode(){
  return &(ynodepool[ypoolptr++]);
}

int phremit(int i){
  //output y's
  ynode_t *yptr;
  ynode_t *lastzy=NULL;

  for (yptr=ylist_head->next; yptr!=ylist_tail; yptr=yptr->next){

    int fval=(fcounts[yptr->uptr->mval] - fcounts[yptr->lptr->mval-1]) - (ecounts[i+1] - ecounts[yptr->yval]);
    cout<<"when i="<<i<<", "<<"y:"<<yptr->yval<<" "<<"fval:"<<fval<<endl;	
    cout<<"fcounts["<<yptr->uptr->mval<<"]-fcounts["<<yptr->lptr->mval-1<<"]"<<endl;
    if (fval==0){
      cout<<"["<<yptr->yval<<", "<<i<<"]"<<endl;

      lparens[yptr->yval]++;
      rparens[i]++;
      lastzy=yptr;
    }
    else
      break;
  }

#ifdef RECURSIVE
  if (lastzy){
    lastzy->prev=ylist_head;
    ylist_head->next=lastzy;
      

    //updating ulist, llist
    ulist_head->next=lastzy->uptr;
    lastzy->uptr->prev=ulist_head;
    
    llist_head->next=lastzy->lptr;
    lastzy->lptr->prev=llist_head;     
  }
#endif
}

int bracketing(){  
  mnode_t *uptr;
  mnode_t *lptr;

  ynode_t *yptr;
  ynode_t *ystarptr;
  ynode_t *yplusptr;

  int prevu, prevl;
  int firstnonempty;

  lparens.resize(n+1);
  lparens.clear();
  rparens.resize(n+1);
  rparens.clear();

  mpoolptr=0;
  ypoolptr=0;

  //initialize, create heads and tails, link them together
  ulist_head=getmnode();
  llist_head=getmnode();
  ylist_head=getynode();

  ulist_head->prev=NULL;
  llist_head->prev=NULL;
  ylist_head->prev=NULL;

  ulist_tail=getmnode();
  llist_tail=getmnode();
  ylist_tail=getynode();

  ulist_tail->next=NULL;
  llist_tail->next=NULL;
  ylist_tail->next=NULL;

  ulist_head->next=getmnode();
  llist_head->next=getmnode();
  ylist_head->next=getynode();


  ulist_head->next->next=ulist_tail;
  ulist_head->next->prev=ulist_head;
  ulist_tail->prev=ulist_head->next;


  llist_head->next->next=llist_tail;
  llist_head->next->prev=llist_head;
  llist_tail->prev=llist_head->next;

  ylist_head->next->next=ylist_tail;
  ylist_head->next->prev=ylist_head;
  ylist_tail->prev=ylist_head->next;

  //cross-linking between u, l, and y
  ulist_head->next->ly=ylist_head->next;
  llist_head->next->ly=ylist_head->next;
  
  ylist_head->next->uptr=ulist_head->next;
  ylist_head->next->lptr=llist_head->next;


  //skip the empty columns
  for (firstnonempty=0; firstnonempty<n && alignment[firstnonempty].size()==0; firstnonempty++);

  if (firstnonempty<n){
    //first column 
    ylist_head->next->yval=firstnonempty;
    ylist_head->next->uptr->mval=u[firstnonempty];
    ylist_head->next->lptr->mval=l[firstnonempty];

    prevu=u[firstnonempty];
    prevl=l[firstnonempty];

    cout<<"u["<<firstnonempty<<"]:"<<prevu<<endl;
    cout<<"l["<<firstnonempty<<"]:"<<prevl<<endl;

    phremit(firstnonempty);
  }


  //going from left to right
  for (int i=firstnonempty+1; i<n; i++) {
    //append a ynode to ylist
     
    //skip the words that are not aligned to anything
    if (alignment[i].size()==0)
      continue;
    
    //otherwise, we have at least one word in each column

    //append a y node to the right end of y list
    ynode *py=getynode();
    py->yval=i;

    py->prev=ylist_head;
    py->next=ylist_head->next;
    ylist_head->next=py;
    py->next->prev=py;

    //link to the old u,l
    py->uptr=ulist_head->next;
    py->lptr=llist_head->next;

    cout<<"u["<<i<<"]:"<<u[i]<<endl;
    cout<<"l["<<i<<"]:"<<l[i]<<endl;

    cout<<"{{{ ";
    for (ynode *itry=ylist_head->next; itry!=ylist_tail; itry=itry->next){
      cout<<itry->yval<<" ";
    }
    cout<<" }}}"<<endl;

    if (u[i]<=prevu){
      //simply append a unode
      mnode *pu=getmnode();
      
      pu->prev=ulist_head;
      pu->next=ulist_head->next;
      ulist_head->next=pu;
      pu->next->prev=pu;

      pu->mval=u[i];
      pu->ly=py;
      
      py->uptr=pu;

    }

    if (l[i]>=prevl){
      //simply append a l node
      mnode *pl=getmnode();
      
      pl->prev=llist_head;
      pl->next=llist_head->next;
      llist_head->next=pl;
      pl->next->prev=pl;

      pl->mval=l[i];
      pl->ly=py;
      
      py->lptr=pl;      
    }
    /////////////////////////
    //upward
    //    if (u[i]>=prevu){
    {
      
      //trace back to y*
      uptr=ulist_head->next;
      cout<<"prevu"<<prevu<<endl;

      if (u[i]>prevu){
	while (uptr->next && uptr->mval<u[i]){
	  uptr=uptr->next;
	}
	
	//found an existing higher step, or from the end node down to the last node
	uptr=uptr->prev; //come down to y* level
      }


      //raise ulist at y*
      uptr->mval=u[i]; //raising
      ylist_head->next->uptr=uptr;
      

      //lemma 4.1
      //remove all the lower unodes

      if (uptr->prev!=ulist_head){//we have some lower steps that can be removed
	//remove the corresponding y's
	ylist_head->next->next=uptr->prev->ly->next;
	uptr->prev->ly->next->prev=ylist_head->next;

	//relink the ulist
	uptr->prev=ulist_head;
	ulist_head->next=uptr;

	if (l[i]>=prevl){// a new lnode has been appended
	  if (llist_head->next!=ylist_head->next->next->lptr){
	    //relink the llist

	    llist_head->next->next=ylist_head->next->next->lptr;
	    ylist_head->next->next->lptr->prev=llist_head->next;
	  }
	}
 	else{ // no new lnode has been appended
 	  llist_head->next=ylist_head->next->next->lptr;
 	  ylist_head->next->next->lptr->prev=llist_head;
 	}
	
	cout<<"lemma 4.1 "<<i<<endl;
      }
      

      //lemma 4.2
      //iteratively remove all yi's that do not satisfy the monotonically increasing property
      ystarptr=uptr->ly;

      cout<<"i="<<i<<","<<"ystar"<<ystarptr->yval<<endl;

      if (ystarptr->next->next){//the next y exists
	yplusptr=ystarptr->next;
	
	bool mono=true;
	do{//function f evaluated here 
	  int fstar=(fcounts[ystarptr->uptr->mval] - fcounts[ystarptr->lptr->mval-1]) - (ecounts[i+1] - ecounts[ystarptr->yval]);
	  int fplus=(fcounts[yplusptr->uptr->mval] - fcounts[yplusptr->lptr->mval-1]) - (ecounts[i+1] - ecounts[yplusptr->yval]);

	  if (fstar > fplus){
	    mono=false;
	    //removing, relinking
	    ystarptr->prev->next=yplusptr;
	    yplusptr->prev=ystarptr->prev;

	    //tracing back
	    ystarptr=ystarptr->prev;

	    cout<<"lemma 4.2 "<<i<<endl;
	  }
	  else
	    mono=true;

	}while (!mono && ystarptr!=ylist_head);

	//updating ulist, llist
	if (ystarptr==ylist_head){//every y in front of y* has been removed. is that possible?
	  cout<<"right end"<<endl;
	  ulist_head->next=yplusptr->uptr;
	  yplusptr->uptr->prev=ulist_head;

	  llist_head->next=yplusptr->lptr;
	  yplusptr->lptr->prev=llist_head; 

	}
	else{

	//relinking ulist and llist
	if (yplusptr->uptr->mval > ystarptr->uptr->mval){
	  ystarptr->uptr->ly=ystarptr;
	  ystarptr->uptr->next=yplusptr->uptr;
	  yplusptr->uptr->prev=ystarptr->uptr;
	}
	if (yplusptr->lptr->mval < ystarptr->lptr->mval){
	  ystarptr->lptr->ly=ystarptr;
	  ystarptr->lptr->next=yplusptr->lptr;
	  yplusptr->lptr->prev=ystarptr->lptr;
	}
	
	}
      }

    }//end of upward case
    ///////////////////////////////////

    if (ystarptr==ylist_head){
      //append a y node to the right end of y list
      py=getynode();
      py->yval=i;

      py->prev=ylist_head;
      py->next=ylist_head->next;
      ylist_head->next=py;
      py->next->prev=py;

      //link to the old u,l
      py->uptr=ulist_head->next;
      py->lptr=llist_head->next;

    }

    //downward
    //    if (l[i]<=prevl){


    //    if (ystarptr!=ylist_head)
    {

      prevl=ylist_head->next->next->lptr->mval;


      if (llist_head->next->mval<l[i]){
	//simply append a l node
	mnode *pl=getmnode();
      
	pl->prev=llist_head;
	pl->next=llist_head->next;
	llist_head->next=pl;
	pl->next->prev=pl;
	
	pl->mval=l[i];
	pl->ly=py;
	  
	py->lptr=pl;      
      }

      //trace back to y* on llist
      lptr=llist_head->next;
      
      if (l[i]<prevl){
	while (lptr->next && lptr->mval>l[i]){
	  lptr=lptr->next;
	}	
	
	//found an existing lower step, or from the end node down to the last node
	lptr=lptr->prev; //come down to y* level

      }

      //lower llist at y*
      lptr->mval=l[i]; //raising

      ylist_head->next->lptr=lptr;


      //lemma 4.1
      //remove all the higher lnodes
	
      if (lptr->prev!=llist_head){//we have some lower steps that can be removed
	//remove the corresponding y's
	cout<<"from "<<lptr->prev->ly->yval<<" to "<<ylist_head->next->next->yval<<endl;
	ylist_head->next->next=lptr->prev->ly->next;
	lptr->prev->ly->next->prev=ylist_head->next;  //skipping the newly inserted yi

	//relink the llist
	lptr->prev=llist_head;
	llist_head->next=lptr;

	//relink the ulist
	if (ulist_head->next!=ylist_head->next->next->uptr){
	  ulist_head->next->next=ylist_head->next->next->uptr;
	  ylist_head->next->next->uptr->prev=ulist_head->next;
	  ulist_head->next->ly=ylist_head->next;	  
	}

	cout<<"lemma 4.1 "<<i<<endl;
      }
      

      //lemma 4.2
      //iteratively remove all yi's that do not satisfy the monotonically increasing property
      ystarptr=lptr->ly;
      cout<<"i="<<i<<","<<"ystar"<<ystarptr->yval<<endl;

      if (ystarptr->next->next){//the next y exists
	yplusptr=ystarptr->next;
	
	bool mono=true;
	do{//evaluation f
	  int fstar=(fcounts[ystarptr->uptr->mval] - fcounts[ystarptr->lptr->mval-1]) - (ecounts[i+1] - ecounts[ystarptr->yval]);
	  int fplus=(fcounts[yplusptr->uptr->mval] - fcounts[yplusptr->lptr->mval-1]) - (ecounts[i+1] - ecounts[yplusptr->yval]);


	  if (fstar > fplus){
	    mono=false;
	    //removing, relinking
	    ystarptr->prev->next=yplusptr;
	    yplusptr->prev=ystarptr->prev;

	    //tracing back
	    ystarptr=ystarptr->prev;
	    cout<<"lemma 4.2 "<<i<<endl;
	  }
	  else
	    mono=true;

	}while (!mono && ystarptr!=ylist_head);

	//updating ulist, llist
	if (ystarptr==ylist_head){//every y in front of y* has been removed. is that possible?
	  cout<<"right end"<<endl;
	  ulist_head->next=yplusptr->uptr;
	  yplusptr->uptr->prev=ulist_head;

	  llist_head->next=yplusptr->lptr;
	  yplusptr->lptr->prev=llist_head; 


	}
	else{

	//relinking
	if (yplusptr->uptr->mval > ystarptr->uptr->mval){
	  ystarptr->uptr->ly=ystarptr;
	  ystarptr->uptr->next=yplusptr->uptr;
	  yplusptr->uptr->prev=ystarptr->uptr;
	}
	if (yplusptr->lptr->mval < ystarptr->lptr->mval){
	  ystarptr->lptr->ly=ystarptr;
	  ystarptr->lptr->next=yplusptr->lptr;
	  yplusptr->lptr->prev=ystarptr->lptr;
	}
	
	}
      }      

    }

    if (ystarptr!=ylist_head){

      //      cout<<"tail"<<endl;
      ystarptr=ylist_head->next;
      yplusptr=ystarptr->next;
      
      bool mono=true;

      do{//evaluation f
	int fstar=(fcounts[ystarptr->uptr->mval] - fcounts[ystarptr->lptr->mval-1]) - (ecounts[i+1] - ecounts[ystarptr->yval]);
	int fplus=(fcounts[yplusptr->uptr->mval] - fcounts[yplusptr->lptr->mval-1]) - (ecounts[i+1] - ecounts[yplusptr->yval]);


	if (fstar > fplus){
	  mono=false;
	  //removing, relinking
	  ystarptr->prev->next=yplusptr;
	  yplusptr->prev=ystarptr->prev;

	  //tracing back
	  ystarptr=ystarptr->prev;
	  cout<<"tail lemma 4.2 "<<i<<endl;
	}
	else
	  mono=true;

      }while (!mono && ystarptr!=ylist_head);


      //updating ulist, llist
      if (ystarptr==ylist_head){//every y in front of y* has been removed. is that possible?

	ulist_head->next=yplusptr->uptr;
	yplusptr->uptr->prev=ulist_head;

	llist_head->next=yplusptr->lptr;
	yplusptr->lptr->prev=llist_head; 
      }
    }

    phremit(i);

    prevu=ulist_head->next->mval;
    prevl=llist_head->next->mval;
  }


  return 0;
}

int output(){
  for (int i=0; i<=n; i++){
    for (int j=0; j<lparens[i]; j++)
      cerr<<" ( ";
    if (i<n)
      cerr<<i+1<<" ";
    for (int j=0; j<rparens[i]; j++)
      cerr<<" ) ";
  }
  cerr<<endl;
  
  return 0;
}

int totalsen=0;

#define MAXLINE 1000
int main(int argc, char *argv[]){
  char linebuf[MAXLINE];

  while (!cin.eof()) {
    cin.getline(linebuf, MAXLINE);
    string linestr(linebuf);
    istringstream linestream(linestr);
    
    //header, n x m
    n=m=0;
    linestream>>n;
    linestream>>m;

    if (n==0){
      cout<<endl;
      continue;
    }

    //initialize the matrix
    alignment.clear();
    alignment.resize(n);
    ecounts.clear();
    ecounts.resize(n+1);
    fcounts.clear();
    fcounts.resize(m+1);
    u.clear();
    u.resize(n);
    l.clear();
    l.resize(n);

    //read in the word links column by column
    for (int i=0; i<n; i++) {
      cin.getline(linebuf, MAXLINE);
      string linestr(linebuf);
      istringstream linestream(linestr);
      
      //initialize u/l
      u[i]=0;
      l[i]=m+1;

      //read the alignment dots for each column
      int j;
      while (linestream>>j){
	alignment[i].push_back(j);
	ecounts[i+1]++;
	fcounts[j]++;
	if (j>u[i])
	  u[i]=j;
	if (j<l[i])
	  l[i]=j;
      }

    }//end of columns

    //accummulate ecounts, fcounts,
    ecounts[0]=0;
    for (int i=1; i<=n; i++)
      ecounts[i]+=ecounts[i-1];
    fcounts[0]=0;
    for (int j=1; j<=m; j++)
      fcounts[j]+=fcounts[j-1];

    //debugging

//     for (int i=0; i<n; i++){
//       cout<<"ecounts["<<i+1<<"]:"<<ecounts[i+1]<<endl;
//     }
//     for (int i=0; i<m; i++){
//       cout<<"fcounts["<<i+1<<"]:"<<fcounts[i+1]<<endl;
//     }
//     for (int i=0; i<n; i++){
//       cout<<"u["<<i<<"]:"<<u[i]<<", "<<"l["<<i<<"]:"<<l[i]<<endl;
//     }

    cout<<"SENPAIR: "<<++totalsen<<endl;
    bracketing(); //factorize it
    output();
  }
}
