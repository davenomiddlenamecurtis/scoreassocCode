extern "C"
{
#include "scoreassoc.h"
};

// lazily make these global to save allocating them and to avoid stack overflow
float weight[MAX_LOCI],missing_score[MAX_LOCI],func_weight[MAX_LOCI],cc_freq[2][MAX_LOCI],cc_count[2][MAX_LOCI],cc_genocount[2][3][MAX_LOCI];
int rarer[MAX_LOCI];
char names[MAX_LOCI][20],comments[MAX_LOCI][MAX_COMMENT_LENGTH],trios_fn[500];
subject **global_sub;