### POST /api/chat,

## Type Definition

### Selection

> 如果 children 不为空，则视为选区组, 选取组的 polygon 和 indexes 可能为空

```ts
interface Selection {
  id: string; // 选区ID（不能重复）
  emb: string; // 视图名称 umap / spatial
  name: string; // 选区名称 （支持自定义）
  polygon: vec2[]; // 选区顶点坐标
  indexes: number[]; // 选区中细胞索引
  children: Selection[]; // 子选区
}
```

### Record

```ts
type RecordType = "message" | "conversation" | "log";

interface Message {
  id: string;
  role: "user" | "assistant" | "system";
  selectionIds: string[]; // 当前选中的选区Id
  content: string;
  ctime: number;
}

interface Log {
  id: string;
  type: string;
  content: string;
  ctime: number;
}

interface Conversation {
  id: string;
  records: Record[];
  ctime: number;
}

interface Record {
  id: string;
  type: RecordType;
  data: Message | Conversation | Log;
}
```

### Request Body

```ts
interface Request {
  selections: Selection[];
  conversation: Conversation;
}
```

### Response Body

```ts
interface Response {}
```

## Example

### Request

```json
{
  "selections": [
    {
      "id": "region1",
      "emb": "umap",
      "name": "ungrouped",
      "polygon": [],
      "indexes": [],
      "children": [
        {
          "id": "region2",
          "emb": "umap",
          "name": "region2",
          "polygon": [
            [0.3, 0.4],
            [0.4, 0.5],
            [0.6, 0.6],
            [0.1, 0.1]
          ],
          "indexes": [5, 4],
          "children": []
        },
        {
          "id": "2",
          "emb": "spatial",
          "name": "region3",
          "polygon": [
            [0.1, 0.1],
            [0.2, 0.1],
            [0.2, 0.2],
            [0.1, 0.2]
          ],
          "indexes": [12, 23],
          "children": []
        }
      ]
    }
  ],
  "conversation": {
    "id": "conversation-1",
    "records": [
      {
        "id": "record-1",
        "type": "message",
        "data": {
          "id": "message-1",
          "role": "user",
          "selectionIds": ["region1"],
          "content": "Analyze the content of this selection",
          "ctime": 1739186045
        }
      },
      {
        "id": "record-2",
        "type": "message",
        "data": {
          "id": "message-2",
          "role": "assistant",
          "selectionIds": [],
          "content": "Okay, that includes xxx",
          "ctime": 1739186048
        }
      },
      {
        "id": "record-3",
        "type": "conversation",
        "data": {
          "id": "conversation-1",
          "records": [
            {
              "id": "record-4",
              "type": "message",
              "data": {
                "id": "message-3",
                "role": "user",
                "selectionIds": [],
                "content": "hi",
                "ctime": 1739186050
              }
            },
            {
              "id": "record-5",
              "type": "message",
              "data": {
                "id": "message-4",
                "role": "assistant",
                "selectionIds": [],
                "content": "Hello! How can I assist you today?",
                "ctime": 1739186051
              }
            }
          ],
          "ctime": 1739186042
        }
      },
      {
        "id": "record-4",
        "type": "message",
        "data": {
          "id": "message-3",
          "role": "user",
          "selectionIds": ["region2"],
          "content": "Reanalyze the content of this selection",
          "ctime": 1739186052
        }
      }
    ],
    "ctime": 1739186042
  }
}
```

### Response

```json

```
