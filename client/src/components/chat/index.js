import React from "react";
import { connect } from "react-redux";
import { Button, H4, TextArea } from "@blueprintjs/core";

import { Observable } from "rxjs";
import * as globals from "../../globals";

const generateId = (() => {
  let id = 0;

  return () => {
    id += 1;
    return `id-${id}`;
  };
})();

const source = new Observable((observer) => {
  const response =
    "Technology has revolutionized the way we live and work in countless ways over the past few decades. From smartphones that fit in our pockets to self-driving cars, innovations continue to push the boundaries of what is possible. However, with these advancements come new challenges, such as privacy concerns and the need for ethical guidelines in AI development. As technology continues to evolve, it's important for society to consider both its benefits and potential drawbacks.";

  let i = 0;

  const write = () => {
    requestAnimationFrame(() => {
      observer.next(response[i]);

      i += 1;

      if (i < response.length) {
        write();
        return;
      }

      observer.complete();
    });
  };

  write();
});

@connect((state) => ({
  session: state.chat.session,
}))
class Chat extends React.Component {
  constructor(props) {
    super(props);

    this.messagesRef = null;

    this.state = {
      input: "",
      output: "",
    };

    this.bindMessagesRef = this.bindMessagesRef.bind(this);
    this.onChange = this.onChange.bind(this);
    this.onKeyDown = this.onKeyDown.bind(this);
    this.archive = this.archive.bind(this);
  }

  componentDidUpdate(prevProps, prevState) {
    const { session: prevSession } = prevProps;
    const { output: prevOutput } = prevState;

    const { session } = this.props;
    const { output } = this.state;

    if (
      session.messages.length !== prevSession.messages.length ||
      output !== prevOutput
    ) {
      if (this.messagesRef) {
        this.messagesRef.scrollTop = this.messagesRef.scrollHeight;
      }
    }
  }

  onChange(e) {
    this.setState((p) => ({
      ...p,
      input: e.target.value,
    }));
  }

  onKeyDown(e) {
    if (e.key === "Enter" && !e.shiftKey) {
      e.stopPropagation();
      e.preventDefault();
      this.sendMessage();
    }
  }

  bindMessagesRef(ref) {
    this.messagesRef = ref;
  }

  async sendMessage() {
    const { dispatch } = this.props;
    const { input, output } = this.state;
    if (input.length === 0 || output.length > 0) {
      return;
    }
    console.log("input: %s", input);

    dispatch({
      type: "chat: new message",
      message: {
        id: generateId(),
        role: "user",
        content: input,
        timestamp: Date.now(),
      },
    });

    this.setState((p) => ({
      ...p,
      input: "",
    }));

    let chunks = "";

    source.subscribe({
      next: (chunk) => {
        chunks += chunk;
        this.setState((p) => ({
          ...p,
          output: chunks,
        }));
      },
      complete: () => {
        dispatch({
          type: "chat: new message",
          message: {
            id: generateId(),
            role: "bio",
            content: chunks,
            timestamp: Date.now(),
          },
        });
        this.setState((p) => ({
          ...p,
          output: "",
        }));
      },
    });
  }

  archive() {
    const { session } = this.props;
    const { messages } = session;

    const text = messages
      .map(
        (msg) =>
          `[${new Date(msg.timestamp).toLocaleString()}] ${msg.role}: ${
            msg.content
          }`
      )
      .join("\n\n");

    console.log(text);

    const blob = new Blob([text], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "chat_archive.txt";
    a.click();
    URL.revokeObjectURL(url);
  }

  render() {
    const { session } = this.props;

    const { input, output } = this.state;

    const { messages } = session;

    return (
      <div
        style={{
          flex: 2,
          borderTop: `0.5px solid ${globals.lightGrey}`,
          overflowY: "inherit",
          padding: globals.leftSidebarSectionPadding,

          display: "flex",
          flexDirection: "column",
          rowGap: "1em",
        }}
      >
        <div
          style={{
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
          }}
        >
          <H4>Chat</H4>
          <Button
            icon="archive"
            disabled={messages.length === 0 || output.length > 0}
            onClick={this.archive}
          />
        </div>
        <div
          ref={this.bindMessagesRef}
          style={{
            flex: 1,
            overflowY: "auto",
            display: "flex",
            flexDirection: "column",
            rowGap: "1em",
            paddingRight: 16,
          }}
        >
          {messages.map((msg) => (
            <div
              key={msg.id}
              style={{ display: "flex", alignItems: "flex-start", gap: "1ch" }}
            >
              <span style={{ width: "10ch" }}>{msg.role}</span>
              <span
                style={{
                  flex: 1,
                  backgroundColor:
                    msg.role === "user" ? "transparent" : "rgba(0,0,0,0.3)",
                }}
              >
                {msg.content}
              </span>
            </div>
          ))}
          <div
            key="output"
            style={{
              display: output ? "flex" : "none",
              alignItems: "flex-start",
              gap: "1ch",
            }}
          >
            <span style={{ width: "10ch" }}>bio</span>
            <span style={{ flex: 1, backgroundColor: "rgba(0,0,0,0.3)" }}>
              {output}
            </span>
          </div>
        </div>
        <TextArea
          autoFocus
          fill
          style={{ resize: "none", maxHeight: 100 }}
          growVertically
          placeholder="Press Enter to send, Shift+Enter for newline"
          value={input}
          onChange={this.onChange}
          onKeyDown={this.onKeyDown}
        />
      </div>
    );
  }
}

export default Chat;
