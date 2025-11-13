
# MySQL Inner Join Through Junction Tables

To join from A → B → C using through-tables, chain your joins together:

```sql
SELECT a.*, b.*, c.*
FROM A a
INNER JOIN A2B ON a.id = A2B.a_id
INNER JOIN B b ON A2B.b_id = b.id
INNER JOIN B2C ON b.id = B2C.b_id
INNER JOIN C c ON B2C.c_id = c.id;
```

```sql
SELECT a.pmcid, b.id, c.id
FROM env_publications a
INNER JOIN env_publication_2_mutation A2B ON a.id = A2B.env_publication_id
INNER JOIN  antibiotic_res_mutations b ON A2B.antibio_res_mutation_id = b.id
INNER JOIN pdb_2_mutation B2C ON b.id = B2C.antibio_res_mutation_id
INNER JOIN pdb_structures c ON B2C.pdb_id = c.id;

```

## Breaking it down:

1. **Start with A** and alias it
2. **Join A2B** on A's primary key
3. **Join B** using the foreign key from A2B
4. **Join B2C** on B's primary key
5. **Join C** using the foreign key from B2C

## Example with realistic column names:

```sql
SELECT a.name, b.title, c.description
FROM users a
INNER JOIN user_roles ON a.user_id = user_roles.user_id
INNER JOIN roles b ON user_roles.role_id = b.role_id
INNER JOIN role_permissions ON b.role_id = role_permissions.role_id
INNER JOIN permissions c ON role_permissions.permission_id = c.permission_id
WHERE a.user_id = 5;
```

## Key points:

- Each junction table uses **two foreign keys** to connect its neighboring tables
- The join order matters—you can't skip steps
- Use `INNER JOIN` (returns only matching rows) or `LEFT JOIN` if you want to include records without matches
- Alias tables for cleaner, more readable queries


# Django ORM Multi-Table Inner Join

The approach depends on your relationship structure. Here are the main scenarios:

## Scenario 1: Using `select_related()` (ForeignKey chain)

```python
# Get all A records with their B and C data
results = A.objects.select_related(
    'b',        # A → B
    'b__c'      # B → C
)

# Or with filtering
results = A.objects.filter(
    b__c__some_field='value'
).select_related('b', 'b__c')
```

## Scenario 2: Using through-tables (ManyToMany)

If you have explicit through models:

```python
# models.py
class A(models.Model):
    name = models.CharField(max_length=100)

class B(models.Model):
    title = models.CharField(max_length=100)
    a_items = models.ManyToManyField(A, through='A2B')

class C(models.Model):
    description = models.CharField(max_length=100)

class A2B(models.Model):
    a = models.ForeignKey(A, on_delete=models.CASCADE)
    b = models.ForeignKey(B, on_delete=models.CASCADE)

class B2C(models.Model):
    b = models.ForeignKey(B, on_delete=models.CASCADE)
    c = models.ForeignKey(C, on_delete=models.CASCADE)
```

**Query:**

```python
# Get all A's with related B's and C's
results = A.objects.prefetch_related(
    'a2b__b__b2c__c'
)

for a in results:
    for a2b in a.a2b_set.all():
        for b2c in a2b.b.b2c_set.all():
            print(a.name, a2b.b.title, b2c.c.description)
```

## Scenario 3: Using `Prefetch` for complex queries

```python
from django.db.models import Prefetch

# More control over querysets
b_prefetch = Prefetch('a2b', queryset=A2B.objects.select_related('b'))
results = A.objects.prefetch_related(b_prefetch)
```

## Key Django ORM concepts:

| SQL | Django |
|-----|--------|
| INNER JOIN (ForeignKey) | `select_related()` |
| INNER JOIN (ManyToMany/reverse) | `prefetch_related()` |
| Filter across relations | `filter(relation__field='value')` |
| Chain relations | Double underscore `__` |

**Pro tip:** Use `prefetch_related()` for junction tables and `select_related()` for direct ForeignKey relationships to minimize queries.

Does this match your Django model structure?