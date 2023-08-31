TList_addObjectAt");

	    newNode->object = newObject;
	    newNode->next = temp;
	    if (prevNode)
		prevNode->next = newNode;
	    return;
	}
	prevNode = temp;
	Pos--;
    }
    if (Pos >= 0)
	HTList_addObject(prevNode, newObject);

    return;
}

/*	Unlink specified object from list.
 *	It does not free memory.
 */
BOOL HTList_unlinkObject(HTList *me, void *oldObject)
{
    HTList *temp = me;
    HTList *prevNode;

    if (temp && oldObject) {
	while (temp->next) {
	    prevNode = temp;
	    temp = temp->next;
	    if (temp->object == oldObject) {
		prevNode->next = temp->next;
		temp->next = NULL;
		temp->object = NULL;
		return YES;	/* Success */
	    }
	}
    }
    return NO;			/* object not found or NULL list */
}

/*	Remove specified object from list.
*/
BOOL HTList_removeObject(HTList *me, void *oldObject)
{
    HTList *temp = me;
    HTList *prevNode;

    if (temp && oldObject) {
	while (temp->next) {
	    prevNode = temp;
	    temp = temp->next;
	    if (temp->object == oldObject) {
		prevNode->next = temp->next;
		FREE(temp);
		return YES;	/* Success */
	    }
	}
    }
    return NO;			/* object not found or NULL list */
}

/*	Remove object at a given position in the list, where 0 is the
 *	object pointed to by the head (returns a pointer to the element
 *	(->object) for the object, and NULL if the list is empty, or
 *	if it doesn't exist - Yuk!).
 */
void *HTList_removeObjectAt(HTList *me, int position)
{
    HTList *temp = me;
    HTList *prevNode;
    int pos = position;
    void *result = NULL;

    if (temp != NULL && pos >= 0) {
	prevNode = temp;
	while ((temp = temp->next) != NULL) {
	    if (pos == 0) {
		prevNode->next = temp->next;
		result = temp->object;
		FREE(temp);
		break;
	    }
	    prevNode = temp;
	    pos--;
	}
    }

    return result;
}

/*	Unlink object from START of list (the Last one inserted
 *	via HTList_linkObject(), and pointed to by the head).
 *	It does not free memory.
 */
void *HTList_unlinkLastObject(HTList *me)
{
    HTList *lastNode;
    void *lastObject;

    if (me && me->next) {
	lastNode = me->next;
	lastObject = lastNode->object;
	me->next = lastNode->next;
	lastNode->next = NULL;
	lastNode->object = NULL;
	return lastObject;

    } else {			/* Empty list */
	return NULL;
    }
}

/*	Remove object from START of list (the Last one inserted
 *	via HTList_addObject(), and pointed to by the head).
 */
void *HTList_removeLastObject(HTList *me)
{
    HTList *lastNode;
    void *lastObject;

    if (me && me->next) {
	lastNode = me->next;
	lastObject = lastNode->object;
	me->next = lastNode->next;
	FREE(lastNode);
	return lastObject;

    } else {			/* Empty list */
	return NULL;
    }
}

/*	Remove object from END of list (the First one inserted
 *	via HTList_addObject(), and furthest from the head).
 */
void *HTList_removeFirstObject(HTList *me)
{
    HTList *temp = me;
    HTList *prevNode;
    void *firstObject;

    if (!temp)
	return NULL;

    prevNode = temp;
    if (temp->next) {
	while (temp->next) {
	    prevNode = temp;
	    temp = temp->next;
	}
	fi