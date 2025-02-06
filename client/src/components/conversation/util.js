import { Observable } from "rxjs";

export const messageOutputSource = new Observable((observer) => {
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
