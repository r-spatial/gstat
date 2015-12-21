#define Q_BUFFER_SIZE 100 /* something more practical */

typedef struct {
	union {
		QTREE_NODE *n;
		DPOINT *p;
	} u;
	int is_node; /* is u the QTREE_NODE (1) or rather the DPOINT (0) ? */
	double dist2; /* squared distance to target location */
} QUEUE_NODE;

typedef struct q_element {
	struct q_element *next;
	QUEUE_NODE el;
} Q_ELEMENT;

typedef struct {
	int length, max_length;
	Q_ELEMENT
		*head,  /* pointer to first element in queue, NULL if empty */
		*empty; /* pointer to empty elements (a stack), NULL if none left */
	int blocks; /* size of memory block */
	Q_ELEMENT **block; /* pointers to malloc'ed memory blocks */
	int (CDECL *cmp)(const QUEUE_NODE *a, const QUEUE_NODE *b);
	/* qsort-able element comparison function */
} QUEUE;

QUEUE *init_queue(QUEUE *q, int (CDECL *cmp)(const QUEUE_NODE *a, const QUEUE_NODE *b));
QUEUE_NODE dequeue(QUEUE *q);
void enqueue(QUEUE *q, QUEUE_NODE *qpt, int n);
void free_queue(QUEUE *q);
