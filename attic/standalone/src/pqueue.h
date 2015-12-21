#ifndef PRIO_Q_H
# define PRIO_Q_H

#ifndef Q_BUFFER_SIZE
# ifdef QUEUE_MAIN
#  define Q_BUFFER_SIZE 5 /* testing purposes mainly */
# else
#  define Q_BUFFER_SIZE 100 /* something more practical */
# endif
#endif

#ifdef SEARCH_H
typedef struct {
	union {
		QTREE_NODE *n;
		DPOINT *p;
	} u;
	int is_node; /* is u the QTREE_NODE (1) or rather the DPOINT (0) ? */
	double dist2; /* squared distance to target location */
} QUEUE_NODE;
# define Q_ELEMENT_WHAT QUEUE_NODE
#else
# define Q_ELEMENT_WHAT double
#endif

typedef struct q_element {
	struct q_element *next;
#ifdef QUEUE_MAIN
	double el;
#else
	QUEUE_NODE el;
#endif
} Q_ELEMENT;

typedef struct {
	int length, max_length;
	Q_ELEMENT
		*head,  /* pointer to first element in queue, NULL if empty */
		*empty; /* pointer to empty elements (a stack), NULL if none left */
	int blocks; /* size of memory block */
	Q_ELEMENT **block; /* pointers to malloc'ed memory blocks */
	int (CDECL *cmp)(const Q_ELEMENT_WHAT *a, const Q_ELEMENT_WHAT *b);
	/* qsort-compatible element comparison function */
} QUEUE;

QUEUE *init_queue(QUEUE *q,
	int (CDECL *cmp)(const Q_ELEMENT_WHAT *a, const Q_ELEMENT_WHAT *b));
Q_ELEMENT_WHAT dequeue(QUEUE *q);
void enqueue(QUEUE *q, Q_ELEMENT_WHAT *qpt, int n);
void free_queue(QUEUE *q);

#endif /* PRIO_Q_H */
