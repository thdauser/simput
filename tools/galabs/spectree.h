extern int searchlev;

struct speclist{
	char* name;
	char* fullname;
	float emin;
	float emax;
	double integral;
	struct speclist *link[2];
	long nentries;
	float *energy;
	float *pflux;
	int level;
};

struct speclist *skew (struct speclist *root);

struct speclist *split (struct speclist *root);

void nil_init();

struct speclist* return_nil();

struct speclist* speclist_newnode(char* new_name, char* fullname, float e_min, float e_max, double integral, long nsteps, float *evalues, float *pfvalues, int level);

struct speclist* insert_spec(struct speclist* root, char* new_name, char* fullname, float e_min, float e_max, double integral, long nsteps, float *evalues, float *pfvalues);

struct speclist* check_if_exists(struct speclist* list, char* name, int* status);

struct speclist* check_fullname_if_exists(struct speclist* list, char* fullname, int* status);